# TCR clumping detection pipeline.
#
# Identifies neighborhoods in TCR space containing more TCRs than expected
# by chance under a null model of independent VDJ rearrangement.  Inspired
# by the ALICE and TCRnet methods.
#
# All heavy computation is delegated to C++ via Rcpp:
#   - rcpp_calc_background_distributions()  (background estimation)
#   - rcpp_poisson_test_loop()              (Poisson significance testing)
#   - rcpp_tcrdist_radius_neighbors()       (neighbor search)
#
# R functions here orchestrate the pipeline, build results, and perform
# single-linkage clustering.


# ---------------------------------------------------------------------------
# 1. setup_tcr_groups  (exported)
# ---------------------------------------------------------------------------

#' Assign alpha and beta chain group indices
#'
#' For each TCR in \code{tcr_df}, builds a composite key from the chain
#' gene usage and CDR3 sequence, then assigns a 0-based group index so
#' that identical chains share the same group.  These groups are used for
#' same-chain masking during neighbor search and Poisson testing.
#'
#' @param tcr_df A \code{data.frame} with at least the columns \code{va},
#'   \code{ja}, \code{cdr3a}, \code{vb}, \code{jb}, \code{cdr3b}.
#'   Optional columns \code{cdr3a_nucseq}, \code{cdr3b_nucseq}, and
#'   \code{subject_id} are included in the key when present.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{agroups}}{Integer vector of length \code{nrow(tcr_df)}.
#'       0-based alpha-chain group indices.}
#'     \item{\code{bgroups}}{Integer vector of length \code{nrow(tcr_df)}.
#'       0-based beta-chain group indices.}
#'   }
#'
#' @examples
#' \donttest{
#' tcr_df <- data.frame(
#'   va = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
#'   ja = c("TRAJ33*01", "TRAJ33*01", "TRAJ49*01"),
#'   cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
#'   vb = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
#'   jb = c("TRBJ2-7*01", "TRBJ2-7*01", "TRBJ1-1*01"),
#'   cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
#'   stringsAsFactors = FALSE
#' )
#' groups <- setup_tcr_groups(tcr_df)
#' groups$agroups  # c(0, 0, 1)  (first two share alpha chain)
#' }
#' @seealso \code{\link{find_clumping}}
#' @export
setup_tcr_groups <- function(tcr_df) {
    required <- c("va", "ja", "cdr3a", "vb", "jb", "cdr3b")
    missing_cols <- setdiff(required, colnames(tcr_df))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "setup_tcr_groups: tcr_df is missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    n <- nrow(tcr_df)
    if (n == 0L) {
        return(list(agroups = integer(0L), bgroups = integer(0L)))
    }

    # Build composite keys with a control-character separator that will not
    # appear in gene names or sequences.
    sep <- "\x01"

    # Include nucseq columns if present (matches Python tuple structure)
    if ("cdr3a_nucseq" %in% colnames(tcr_df)) {
        alpha_keys <- paste(tcr_df$va, tcr_df$ja, tcr_df$cdr3a,
                            tcr_df$cdr3a_nucseq, sep = sep)
    } else {
        alpha_keys <- paste(tcr_df$va, tcr_df$ja, tcr_df$cdr3a, sep = sep)
    }

    if ("cdr3b_nucseq" %in% colnames(tcr_df)) {
        beta_keys <- paste(tcr_df$vb, tcr_df$jb, tcr_df$cdr3b,
                           tcr_df$cdr3b_nucseq, sep = sep)
    } else {
        beta_keys <- paste(tcr_df$vb, tcr_df$jb, tcr_df$cdr3b, sep = sep)
    }

    # Include subject_id if present
    if ("subject_id" %in% colnames(tcr_df)) {
        alpha_keys <- paste(alpha_keys, tcr_df$subject_id, sep = sep)
        beta_keys  <- paste(beta_keys, tcr_df$subject_id, sep = sep)
    }

    # Sort unique keys and map to 0-based indices
    unique_alpha <- sort(unique(alpha_keys))
    unique_beta  <- sort(unique(beta_keys))

    alpha_map <- stats::setNames(seq_along(unique_alpha) - 1L, unique_alpha)
    beta_map  <- stats::setNames(seq_along(unique_beta) - 1L, unique_beta)

    agroups <- unname(alpha_map[alpha_keys])
    bgroups <- unname(beta_map[beta_keys])

    list(agroups = agroups, bgroups = bgroups)
}


# ---------------------------------------------------------------------------
# 2. .estimate_background_tcrdist_distributions  (internal)
# ---------------------------------------------------------------------------

#' Estimate background paired TCRdist distributions
#'
#' For each foreground TCR, estimates the probability of seeing a paired
#' TCRdist score at or below each integer distance from 0 to \code{max_dist}
#' under a null model where alpha and beta chains are independently drawn
#' from shuffled backgrounds.
#'
#' @param organism Character string. Organism key (e.g. \code{"human"}).
#' @param tcr_df A data.frame with columns \code{va, ja, cdr3a, cdr3a_nucseq,
#'   vb, jb, cdr3b, cdr3b_nucseq}.
#' @param max_dist Integer. Maximum paired TCRdist to consider.
#' @param num_random_samples Integer. Number of random background chains to
#'   generate per chain type.
#' @param pseudocount Numeric. Pseudocount added before normalizing.
#' @param preserve_vj_pairings Logical. Preserve V-J pairings in resampling.
#' @param bg_tcrs Optional data.frame used for background generation.
#' @return Numeric matrix of dimensions \code{nrow(tcr_df) x (max_dist + 1)}.
#' @keywords internal
.estimate_background_tcrdist_distributions <- function(
    organism,
    tcr_df,
    max_dist,
    num_random_samples = 50000L,
    pseudocount = 0.25,
    preserve_vj_pairings = FALSE,
    bg_tcrs = NULL
) {
    max_dist <- as.integer(max_dist + 0.1)

    # ---- Validate nucseq columns -------------------------------------------
    required_nucseq <- c("cdr3a_nucseq", "cdr3b_nucseq")
    missing_nuc <- setdiff(required_nucseq, colnames(tcr_df))
    if (length(missing_nuc) > 0L) {
        stop(sprintf(
            paste0(".estimate_background_tcrdist_distributions: ",
                   "tcr_df is missing required nucleotide sequence columns: %s. ",
                   "These are needed for junction parsing and background generation."),
            paste(missing_nuc, collapse = ", ")
        ))
    }

    # ---- Determine background TCRs ----------------------------------------
    if (is.null(bg_tcrs)) {
        bg_tcrs <- tcr_df
    }

    # Optimize V/J alleles for better nucleotide matching
    tcrs_for_bg <- .find_alternate_alleles_for_tcrs(organism, bg_tcrs,
                                                     verbose = FALSE)

    # ---- Parse junctions and resample background chains --------------------
    message("estimate_background_tcrdist_distributions: parsing junctions...")
    junctions_df <- .parse_tcr_junctions(organism, tcrs_for_bg)

    message("estimate_background_tcrdist_distributions: resampling ",
            num_random_samples, " background alpha chains...")
    bg_alpha <- .resample_shuffled_tcr_chains(
        organism, num_random_samples, "A", junctions_df,
        preserve_vj_pairings = preserve_vj_pairings)

    message("estimate_background_tcrdist_distributions: resampling ",
            num_random_samples, " background beta chains...")
    bg_beta <- .resample_shuffled_tcr_chains(
        organism, num_random_samples, "B", junctions_df,
        preserve_vj_pairings = preserve_vj_pairings)

    if (is.null(bg_alpha) || nrow(bg_alpha) == 0L) {
        stop(".estimate_background_tcrdist_distributions: ",
             "failed to generate any background alpha chains")
    }
    if (is.null(bg_beta) || nrow(bg_beta) == 0L) {
        stop(".estimate_background_tcrdist_distributions: ",
             "failed to generate any background beta chains")
    }

    # ---- Compute background distributions via C++ --------------------------
    message("estimate_background_tcrdist_distributions: computing ",
            "background distributions for ", nrow(tcr_df), " TCRs...")

    v_dist_a <- .compute_v_region_distance_matrix(organism, "A")
    v_dist_b <- .compute_v_region_distance_matrix(organism, "B")

    tcrdist_freqs <- rcpp_calc_background_distributions(
        fg_va    = as.character(tcr_df$va),
        fg_cdr3a = as.character(tcr_df$cdr3a),
        fg_vb    = as.character(tcr_df$vb),
        fg_cdr3b = as.character(tcr_df$cdr3b),
        bg_va    = as.character(bg_alpha$v_gene),
        bg_cdr3a = as.character(bg_alpha$cdr3),
        bg_vb    = as.character(bg_beta$v_gene),
        bg_cdr3b = as.character(bg_beta$cdr3),
        v_dist_a = v_dist_a,
        v_dist_b = v_dist_b,
        max_dist = max_dist,
        pseudocount = pseudocount
    )

    message("estimate_background_tcrdist_distributions: done.")
    tcrdist_freqs
}


# ---------------------------------------------------------------------------
# 3. .single_linkage_clumping  (internal)
# ---------------------------------------------------------------------------

#' Single-linkage clustering of clumped TCRs
#'
#' Given a set of clumped TCR indices and their neighbor sets, performs
#' single-linkage clustering by iteratively propagating the smallest
#' neighbor index until convergence.
#'
#' @param all_clumped_nbrs Named list (character keys are 0-based clone
#'   indices). Each value is an integer vector of 0-based neighbor indices.
#' @param num_clones Integer. Total number of clones.
#' @param is_clumped Logical vector of length \code{num_clones}.
#' @return Integer vector of length \code{num_clones}. 0 means not
#'   clumped; positive integers are clumping group IDs (reordered by
#'   decreasing group size).
#' @keywords internal
.single_linkage_clumping <- function(all_clumped_nbrs, num_clones, is_clumped) {
    clumped_inds <- sort(as.integer(names(all_clumped_nbrs)))

    # Make neighbor sets symmetric
    for (ii_str in names(all_clumped_nbrs)) {
        ii <- as.integer(ii_str)
        for (nbr in all_clumped_nbrs[[ii_str]]) {
            nbr_str <- as.character(nbr)
            if (!is.null(all_clumped_nbrs[[nbr_str]])) {
                all_clumped_nbrs[[nbr_str]] <- unique(c(
                    all_clumped_nbrs[[nbr_str]], ii))
            }
        }
    }

    # Initialize: each clumped index points to its smallest neighbor
    all_smallest_nbr <- list()
    for (ii_str in names(all_clumped_nbrs)) {
        all_smallest_nbr[[ii_str]] <- min(all_clumped_nbrs[[ii_str]])
    }

    # Iterate until convergence
    repeat {
        updated <- FALSE
        for (ii_str in names(all_smallest_nbr)) {
            current_nbr <- all_smallest_nbr[[ii_str]]
            nbrs <- all_clumped_nbrs[[ii_str]]
            candidates <- vapply(as.character(nbrs), function(x) {
                all_smallest_nbr[[x]]
            }, integer(1L))
            new_nbr <- min(current_nbr, min(candidates))
            if (new_nbr != current_nbr) {
                all_smallest_nbr[[ii_str]] <- new_nbr
                updated <- TRUE
            }
        }
        if (!updated) break
    }

    # Assign cluster numbers
    clusters <- integer(num_clones)  # 0 = not clumped

    cluster_number <- 0L
    cluster_sizes <- list()

    for (ii in clumped_inds) {
        ii_str <- as.character(ii)
        nbr <- all_smallest_nbr[[ii_str]]
        if (ii == nbr) {
            cluster_number <- cluster_number + 1L
            members <- vapply(names(all_smallest_nbr), function(x) {
                all_smallest_nbr[[x]] == nbr
            }, logical(1L))
            member_ids <- as.integer(names(all_smallest_nbr)[members])
            # 1-based indexing for 0-based clone IDs
            clusters[member_ids + 1L] <- cluster_number
            cluster_sizes[[as.character(cluster_number)]] <- length(member_ids)
        }
    }

    # Reorder clusters by decreasing size
    if (cluster_number > 0L) {
        sizes <- unlist(cluster_sizes)
        size_order <- order(sizes, decreasing = TRUE)
        remap <- integer(cluster_number)
        for (i in seq_along(size_order)) {
            remap[size_order[i]] <- i
        }
        remap <- c(0L, remap)  # index 1 maps cluster 0 (not clumped) -> 0

        new_clusters <- integer(num_clones)
        for (idx in seq_len(num_clones)) {
            old_id <- clusters[idx]
            new_clusters[idx] <- remap[old_id + 1L]
        }
        clusters <- new_clusters
    }

    clusters
}


# ---------------------------------------------------------------------------
# 4. .empty_clumping_results  (internal helper)
# ---------------------------------------------------------------------------

#' Create an empty clumping results data.frame
#' @return A zero-row data.frame with the correct columns.
#' @keywords internal
.empty_clumping_results <- function() {
    data.frame(
        clump_type          = character(0L),
        clone_index         = integer(0L),
        nbr_radius          = integer(0L),
        pvalue_adj          = numeric(0L),
        num_nbrs            = integer(0L),
        expected_num_nbrs   = numeric(0L),
        raw_count           = numeric(0L),
        va                  = character(0L),
        ja                  = character(0L),
        cdr3a               = character(0L),
        vb                  = character(0L),
        jb                  = character(0L),
        cdr3b               = character(0L),
        clumping_group      = integer(0L),
        clonotype_fdr_value = numeric(0L),
        stringsAsFactors    = FALSE
    )
}


# ---------------------------------------------------------------------------
# 5. find_clumping  (exported)
# ---------------------------------------------------------------------------

#' Find TCR clumping in TCR space
#'
#' Identifies neighborhoods in TCR space containing more TCRs than expected
#' by chance under a null model of independent VDJ rearrangement.  For each
#' TCR, counts how many other TCRs fall within a set of fixed TCRdist radii,
#' and compares the observed count to the Poisson expectation derived from
#' background distributions.
#'
#' The pipeline:
#' \enumerate{
#'   \item Estimate per-TCR background frequency distributions via shuffled
#'     chain resampling.
#'   \item Assign alpha/beta chain groups for same-chain masking.
#'   \item Find all neighbors within \code{max(radii)} using
#'     \code{\link{tcrdist_radius_neighbors}}.
#'   \item Run Poisson tests at each radius (C++ via
#'     \code{rcpp_poisson_test_loop}).
#'   \item Perform single-linkage clustering of significant clumps.
#' }
#'
#' @param tcr_df A \code{data.frame} with columns \code{va}, \code{ja},
#'   \code{cdr3a}, \code{cdr3a_nucseq}, \code{vb}, \code{jb}, \code{cdr3b},
#'   \code{cdr3b_nucseq}.
#' @param organism Character string. Organism key (e.g. \code{"human"},
#'   \code{"mouse"}).
#' @param radii Integer vector. TCRdist radii to test.
#'   Default \code{c(24L, 48L, 72L, 96L)}.
#' @param num_random_samples Integer. Number of random background chains
#'   per chain type. Default \code{50000L}.
#' @param pvalue_threshold Numeric. Maximum adjusted p-value to include
#'   in results. Default \code{1.0} (include all).
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#' @param clusters_gex Integer vector of length \code{nrow(tcr_df)}, or
#'   \code{NULL}. If provided, also tests for TCR clumps within each GEX
#'   cluster. Default \code{NULL}.
#' @param bg_tcrs Optional data.frame. If provided, used for background
#'   generation instead of \code{tcr_df}. Default \code{NULL}.
#' @param preserve_vj_pairings Logical. Preserve V-J pairings in background
#'   resampling. Default \code{FALSE}.
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{\code{results_df}}{A data.frame sorted by \code{pvalue_adj} with
#'       columns: \code{clump_type}, \code{clone_index} (0-based),
#'       \code{nbr_radius}, \code{pvalue_adj}, \code{num_nbrs},
#'       \code{expected_num_nbrs}, \code{raw_count}, \code{va}, \code{ja},
#'       \code{cdr3a}, \code{vb}, \code{jb}, \code{cdr3b},
#'       \code{clumping_group}, \code{clonotype_fdr_value}.}
#'     \item{\code{is_clumped}}{Logical vector of length \code{nrow(tcr_df)}.}
#'     \item{\code{clusters}}{Integer vector of length \code{nrow(tcr_df)}.
#'       0 = not clumped, positive = cluster ID.}
#'     \item{\code{all_raw_pvalues}}{Numeric matrix
#'       (\code{nrow(tcr_df)} x \code{length(radii)}).}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- find_clumping(tcr_df, "human")
#' result$results_df
#' sum(result$is_clumped)
#' }
#' @seealso \code{\link{setup_tcr_groups}}, \code{\link{tcrdist_radius_neighbors}}, \code{\link{find_meta_clonotypes}}
#' @importFrom stats ppois p.adjust
#' @export
find_clumping <- function(
    tcr_df,
    organism,
    radii = c(24L, 48L, 72L, 96L),
    num_random_samples = 50000L,
    pvalue_threshold = 1.0,
    verbose = TRUE,
    clusters_gex = NULL,
    bg_tcrs = NULL,
    preserve_vj_pairings = FALSE
) {
    # ---- Input validation --------------------------------------------------
    if (!is.data.frame(tcr_df)) {
        tcr_df <- as.data.frame(tcr_df, stringsAsFactors = FALSE)
    }

    required_cols <- c("va", "ja", "cdr3a", "cdr3a_nucseq",
                       "vb", "jb", "cdr3b", "cdr3b_nucseq")
    missing_cols <- setdiff(required_cols, colnames(tcr_df))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "find_clumping: tcr_df is missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    num_clones <- nrow(tcr_df)
    if (num_clones == 0L) {
        return(list(
            results_df     = .empty_clumping_results(),
            is_clumped     = logical(0L),
            clusters       = integer(0L),
            all_raw_pvalues = matrix(numeric(0L), nrow = 0L, ncol = length(radii))
        ))
    }

    radii <- as.integer(radii + 0.1)
    num_radii <- length(radii)

    if (!is.null(clusters_gex)) {
        if (length(clusters_gex) != num_clones) {
            stop(sprintf(
                "find_clumping: clusters_gex length (%d) != nrow(tcr_df) (%d)",
                length(clusters_gex), num_clones
            ))
        }
        clusters_gex <- as.integer(clusters_gex)
    }

    # ---- Step 1: Estimate background distributions -------------------------
    if (verbose) {
        message("find_clumping: estimating background distributions...")
    }
    bg_freqs <- .estimate_background_tcrdist_distributions(
        organism = organism,
        tcr_df = tcr_df,
        max_dist = max(radii),
        num_random_samples = num_random_samples,
        preserve_vj_pairings = preserve_vj_pairings,
        bg_tcrs = bg_tcrs
    )

    # ---- Step 2: Get alpha/beta chain groups -------------------------------
    groups <- setup_tcr_groups(tcr_df)
    agroups <- groups$agroups  # 0-based
    bgroups <- groups$bgroups  # 0-based

    # ---- Step 3: Find neighbors within max(radii) --------------------------
    if (verbose) {
        message("find_clumping: finding neighbors within radius ",
                max(radii), "...")
    }
    # tcrdist_radius_neighbors returns 1-based indices in $indices
    all_nbr_raw <- tcrdist_radius_neighbors(
        tcr_df, organism, radius = max(radii),
        agroups = agroups, bgroups = bgroups
    )

    # ---- Step 4: Convert neighbor data for C++ (0-based indices) -----------
    # all_nbr_raw[[i]]$indices is 1-based; C++ expects 0-based
    all_nbr_indices <- lapply(all_nbr_raw, function(x) {
        as.integer(x$indices - 1L)
    })
    all_nbr_distances <- lapply(all_nbr_raw, function(x) {
        as.numeric(x$distances)
    })

    # ---- Step 5: Poisson test for each TCR and radius ----------------------
    n_bg_pairs <- as.double(num_random_samples) * as.double(num_random_samples)

    poisson_result <- rcpp_poisson_test_loop(
        all_nbr_indices   = all_nbr_indices,
        all_nbr_distances = all_nbr_distances,
        bg_freqs          = bg_freqs,
        agroups           = agroups,
        bgroups           = bgroups,
        radii             = radii,
        n_bg_pairs        = n_bg_pairs,
        pvalue_threshold  = pvalue_threshold,
        num_clones        = num_clones,
        clusters_gex_nullable = clusters_gex,
        use_conservative_pvalues = TRUE
    )

    is_clumped <- poisson_result$is_clumped
    all_raw_pvalues <- poisson_result$all_raw_pvalues

    # ---- Step 6: Build results data.frame ----------------------------------
    if (length(poisson_result$clone_index) == 0L) {
        if (verbose) {
            message("find_clumping: no significant clumps found.")
        }
        return(list(
            results_df      = .empty_clumping_results(),
            is_clumped      = is_clumped,
            clusters        = integer(num_clones),
            all_raw_pvalues = all_raw_pvalues
        ))
    }

    # clone_index from C++ is 0-based; use +1 for R indexing into tcr_df
    idx_1 <- poisson_result$clone_index + 1L

    results_df <- data.frame(
        clump_type        = poisson_result$clump_type,
        clone_index       = poisson_result$clone_index,
        nbr_radius        = poisson_result$nbr_radius,
        pvalue_adj        = poisson_result$pvalue_adj,
        num_nbrs          = poisson_result$num_nbrs,
        expected_num_nbrs = poisson_result$expected_num_nbrs,
        raw_count         = poisson_result$raw_count,
        va                = tcr_df$va[idx_1],
        ja                = tcr_df$ja[idx_1],
        cdr3a             = tcr_df$cdr3a[idx_1],
        vb                = tcr_df$vb[idx_1],
        jb                = tcr_df$jb[idx_1],
        cdr3b             = tcr_df$cdr3b[idx_1],
        stringsAsFactors  = FALSE
    )

    if (verbose) {
        for (r in seq_len(nrow(results_df))) {
            row <- results_df[r, ]
            prefix <- if (row$clump_type == "global") {
                "tcr_nbrs_global"
            } else {
                "tcr_nbrs_intra"
            }
            message(sprintf(
                "%s: %2d %9.6f radius: %2d pval: %9.1e %9.1f tcr: %s %s %s %s %s %s",
                prefix, row$num_nbrs, row$expected_num_nbrs,
                row$nbr_radius, row$pvalue_adj, row$raw_count,
                row$va, row$ja, row$cdr3a, row$vb, row$jb, row$cdr3b
            ))
        }
    }

    # ---- Step 7: Compute FDR values (Benjamini-Hochberg) -------------------
    fdr_values_flat <- p.adjust(as.vector(all_raw_pvalues), method = "BH")
    fdr_values_mat <- matrix(fdr_values_flat, nrow = num_clones, ncol = num_radii)
    fdr_per_clone <- apply(fdr_values_mat, 1L, min)

    results_df$clonotype_fdr_value <- ifelse(
        results_df$clump_type == "global",
        fdr_per_clone[results_df$clone_index + 1L],
        NA_real_
    )

    # ---- Step 8: Build clumped neighbor sets for clustering ----------------
    all_clumped_nbrs <- list()

    for (row_idx in seq_len(nrow(results_df))) {
        ii_0 <- results_df$clone_index[row_idx]
        ii <- ii_0 + 1L  # 1-based for indexing into all_nbr_* lists
        radius <- results_df$nbr_radius[row_idx]

        # all_nbr_indices are already 0-based
        ii_nbr_indices   <- all_nbr_indices[[ii]]
        ii_nbr_distances <- all_nbr_distances[[ii]]

        # Neighbors within this row's radius that are also clumped
        within_mask <- ii_nbr_distances <= radius
        candidate_nbrs <- ii_nbr_indices[within_mask]
        clumped_nbrs <- candidate_nbrs[is_clumped[candidate_nbrs + 1L]]
        clumped_nbrs <- unique(c(clumped_nbrs, ii_0))

        ii_str <- as.character(ii_0)
        if (!is.null(all_clumped_nbrs[[ii_str]])) {
            all_clumped_nbrs[[ii_str]] <- unique(c(
                all_clumped_nbrs[[ii_str]], clumped_nbrs))
        } else {
            all_clumped_nbrs[[ii_str]] <- clumped_nbrs
        }
    }

    # ---- Step 9: Single-linkage clustering ---------------------------------
    clusters <- .single_linkage_clumping(all_clumped_nbrs, num_clones,
                                          is_clumped)

    results_df$clumping_group <- clusters[results_df$clone_index + 1L]

    # ---- Sort by p-value and return ----------------------------------------
    results_df <- results_df[order(results_df$pvalue_adj), , drop = FALSE]
    rownames(results_df) <- NULL

    list(
        results_df      = results_df,
        is_clumped      = is_clumped,
        clusters        = clusters,
        all_raw_pvalues = all_raw_pvalues
    )
}
