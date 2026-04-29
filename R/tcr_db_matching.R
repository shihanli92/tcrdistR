# TCR database matching: distance-based matching of query TCRs against
# literature TCR databases with background-corrected p-values.
#
# Ported from Python CoNGA: conga/tcr_clumping.py lines 551-920


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Normalize TCR database column names
#'
#' Renames \code{va_gene} to \code{va}, \code{vb_gene} to \code{vb}, etc.
#'
#' @param df A data.frame.
#' @return The same data.frame with normalized column names.
#' @keywords internal
.normalize_tcr_db_columns <- function(df) {
    nms <- colnames(df)
    if (!"va" %in% nms && "va_gene" %in% nms) {
        colnames(df)[colnames(df) == "va_gene"] <- "va"
    }
    if (!"vb" %in% nms && "vb_gene" %in% nms) {
        colnames(df)[colnames(df) == "vb_gene"] <- "vb"
    }
    df
}


#' Empty result data.frame for TCR database matching
#' @keywords internal
.empty_tcr_db_match_df <- function() {
    data.frame(
        tcrdist = integer(0), pvalue_adj = numeric(0),
        fdr_value = numeric(0), query_index = integer(0),
        db_index = integer(0), va = character(0),
        cdr3a = character(0), vb = character(0),
        cdr3b = character(0), stringsAsFactors = FALSE
    )
}


# ---------------------------------------------------------------------------
# find_significant_tcrdist_matches
# ---------------------------------------------------------------------------

#' Find significant TCRdist matches between query and database TCRs
#'
#' Computes paired TCRdist distances between query and database TCRs and
#' converts them to p-values adjusted for both the number of query and
#' database TCRs. Background distributions are estimated via the V(D)J
#' rearrangement model in
#' \code{.estimate_background_tcrdist_distributions()}.
#'
#' @param query_tcrs_df A data.frame with at least columns \code{va} (or
#'   \code{va_gene}), \code{cdr3a}, \code{vb} (or \code{vb_gene}),
#'   \code{cdr3b}. If \code{background_tcrs_df} is NULL, also needs
#'   \code{ja}, \code{jb}, \code{cdr3a_nucseq}, \code{cdr3b_nucseq}.
#' @param db_tcrs_df A data.frame with at least columns \code{va} (or
#'   \code{va_gene}), \code{cdr3a}, \code{vb} (or \code{vb_gene}),
#'   \code{cdr3b}.
#' @param organism Character string (e.g. \code{"human"}, \code{"mouse"}).
#' @param adjusted_pvalue_threshold Numeric. Maximum adjusted p-value to
#'   report. Default \code{1.0}.
#' @param background_tcrs_df Optional data.frame for background generation.
#'   If NULL, uses \code{query_tcrs_df} (which must then include \code{ja},
#'   \code{jb}, and nucseq columns).
#' @param num_random_samples Integer. Number of random background samples.
#'   Default \code{50000L}.
#' @param fixup_alleles Logical. If TRUE, optimize allele assignments in
#'   background TCRs. Default \code{TRUE}.
#' @return A data.frame with columns \code{tcrdist}, \code{pvalue_adj},
#'   \code{fdr_value}, \code{query_index} (0-based), \code{db_index}
#'   (0-based), plus query and db TCR information. Sorted by
#'   \code{pvalue_adj}.
#' @export
find_significant_tcrdist_matches <- function(
    query_tcrs_df, db_tcrs_df, organism,
    adjusted_pvalue_threshold = 1.0,
    background_tcrs_df = NULL,
    num_random_samples = 50000L,
    fixup_alleles = TRUE
) {
    # ---- Normalize column names ---------------------------------------------
    query_tcrs_df <- .normalize_tcr_db_columns(query_tcrs_df)
    db_tcrs_df <- .normalize_tcr_db_columns(db_tcrs_df)

    # ---- Validate required columns ------------------------------------------
    for (col in c("va", "cdr3a", "vb", "cdr3b")) {
        if (!col %in% colnames(query_tcrs_df)) {
            stop("find_significant_tcrdist_matches: query_tcrs_df missing ",
                 "column: ", col)
        }
        if (!col %in% colnames(db_tcrs_df)) {
            stop("find_significant_tcrdist_matches: db_tcrs_df missing ",
                 "column: ", col)
        }
    }

    # ---- Default background -------------------------------------------------
    if (is.null(background_tcrs_df)) {
        background_tcrs_df <- query_tcrs_df
    } else {
        background_tcrs_df <- .normalize_tcr_db_columns(background_tcrs_df)
    }

    bg_required <- c("va", "ja", "cdr3a", "cdr3a_nucseq",
                      "vb", "jb", "cdr3b", "cdr3b_nucseq")
    bg_missing <- setdiff(bg_required, colnames(background_tcrs_df))
    if (length(bg_missing) > 0L) {
        stop("find_significant_tcrdist_matches: background_tcrs_df missing ",
             "columns: ", paste(bg_missing, collapse = ", "))
    }

    nq <- nrow(query_tcrs_df)
    ndb <- nrow(db_tcrs_df)
    num_comparisons <- as.double(nq) * as.double(ndb)

    message("find_significant_tcrdist_matches: num_comparisons: ",
            num_comparisons, " (", nq, " x ", ndb, ")")

    # ---- Background distributions -------------------------------------------
    max_dist <- 200L

    bg_freqs <- .estimate_background_tcrdist_distributions(
        organism = organism,
        tcr_df = query_tcrs_df,
        max_dist = max_dist,
        num_random_samples = num_random_samples,
        bg_tcrs = background_tcrs_df
    )

    # ---- Determine max distance worth matching ------------------------------
    adjusted_bg_freqs <- num_comparisons * bg_freqs
    could_match <- apply(adjusted_bg_freqs, 2,
                          function(col) any(col <= adjusted_pvalue_threshold))

    max_dist_for_matching <- 0L
    while (max_dist_for_matching < max_dist &&
           could_match[max_dist_for_matching + 1L]) {
        max_dist_for_matching <- max_dist_for_matching + 1L
    }

    message("find_significant_tcrdist_matches: max_dist_for_matching: ",
            max_dist_for_matching)

    if (max_dist_for_matching == 0L && !could_match[1L]) {
        message("find_significant_tcrdist_matches: no matches possible")
        return(.empty_tcr_db_match_df())
    }

    # ---- Compute rectangular distance matrix (query x db) -------------------
    dist_block <- tcrdist_rect(query_tcrs_df, db_tcrs_df, organism)

    # ---- Find matches within max_dist_for_matching --------------------------
    match_pairs <- which(dist_block <= max_dist_for_matching, arr.ind = TRUE)
    if (nrow(match_pairs) == 0L) {
        return(.empty_tcr_db_match_df())
    }

    match_qi <- match_pairs[, 1L]    # 1-based query index
    match_di <- match_pairs[, 2L]    # 1-based db index
    match_dists <- dist_block[match_pairs]

    # ---- Compute raw p-values and adjusted p-values -------------------------
    raw_pvalues <- numeric(length(match_qi))
    pvalues_adj <- numeric(length(match_qi))
    for (k in seq_along(match_qi)) {
        raw_pvalues[k] <- bg_freqs[match_qi[k], match_dists[k] + 1L]
        pvalues_adj[k] <- num_comparisons * raw_pvalues[k]
    }

    # ---- FDR computation (matching Python exactly) --------------------------
    n_matches <- length(raw_pvalues)
    argsort <- order(raw_pvalues)
    n_pad <- as.integer(num_comparisons) - n_matches
    if (n_pad < 0L) n_pad <- 0L
    sorted_raw <- c(raw_pvalues[argsort], rep(1.0, n_pad))
    fdr_all <- stats::p.adjust(sorted_raw, method = "BH")

    inv_argsort <- integer(n_matches)
    inv_argsort[argsort] <- seq_len(n_matches)
    fdr_values <- fdr_all[inv_argsort]

    # ---- Filter by adjusted_pvalue_threshold and build results ---------------
    keep <- pvalues_adj <= adjusted_pvalue_threshold
    if (!any(keep)) {
        return(.empty_tcr_db_match_df())
    }

    match_qi <- match_qi[keep]
    match_di <- match_di[keep]
    match_dists <- match_dists[keep]
    pvalues_adj <- pvalues_adj[keep]
    fdr_values <- fdr_values[keep]

    result_rows <- vector("list", length(match_qi))
    for (k in seq_along(match_qi)) {
        qi <- match_qi[k]
        di <- match_di[k]
        row <- list(
            tcrdist = match_dists[k],
            pvalue_adj = pvalues_adj[k],
            fdr_value = fdr_values[k],
            query_index = qi - 1L,  # 0-based
            db_index = di - 1L,     # 0-based
            va = query_tcrs_df$va[qi],
            cdr3a = query_tcrs_df$cdr3a[qi],
            vb = query_tcrs_df$vb[qi],
            cdr3b = query_tcrs_df$cdr3b[qi]
        )
        if ("ja" %in% colnames(query_tcrs_df)) {
            row$ja <- query_tcrs_df$ja[qi]
        }
        if ("jb" %in% colnames(query_tcrs_df)) {
            row$jb <- query_tcrs_df$jb[qi]
        }
        for (tag in colnames(db_tcrs_df)) {
            row[[paste0("db_", tag)]] <- db_tcrs_df[[tag]][di]
        }
        result_rows[[k]] <- as.data.frame(row, stringsAsFactors = FALSE)
    }

    results_df <- do.call(rbind, result_rows)
    results_df <- results_df[order(results_df$pvalue_adj), ]
    rownames(results_df) <- NULL
    results_df
}


# ---------------------------------------------------------------------------
# match_tcrs_to_db
# ---------------------------------------------------------------------------

#' Match TCRs to a literature database
#'
#' Wrapper around \code{\link{find_significant_tcrdist_matches}} that loads
#' a default or custom database file.
#'
#' @param tcr_df A data.frame with columns \code{va, ja, cdr3a,
#'   cdr3a_nucseq, vb, jb, cdr3b, cdr3b_nucseq}.
#' @param organism Character string.
#' @param db_tcrs_tsvfile Path to a TSV file with database TCRs. If NULL,
#'   uses the built-in \code{new_paired_tcr_db_for_matching_nr.tsv} (human
#'   only).
#' @param adjusted_pvalue_threshold Numeric. Default \code{1.0}.
#' @param background_tcrs_df Optional background data.frame.
#' @param num_random_samples Integer. Default \code{50000L}.
#' @return A data.frame of significant matches (see
#'   \code{\link{find_significant_tcrdist_matches}}).
#' @export
match_tcrs_to_db <- function(tcr_df, organism,
                              db_tcrs_tsvfile = NULL,
                              adjusted_pvalue_threshold = 1.0,
                              background_tcrs_df = NULL,
                              num_random_samples = 50000L) {
    if (is.null(db_tcrs_tsvfile)) {
        if (organism != "human") {
            warning("match_tcrs_to_db: default paired database only ",
                    "available for human; returning empty data.frame")
            return(.empty_tcr_db_match_df())
        }
        db_tcrs_tsvfile <- system.file(
            "extdata", "new_paired_tcr_db_for_matching_nr.tsv",
            package = "tcrdistR")
        if (!nzchar(db_tcrs_tsvfile)) {
            stop("match_tcrs_to_db: paired TCR database file not found. ",
                 "Ensure the package is properly installed.")
        }
        message("match_tcrs_to_db: matching to default literature TCR database")
    }

    message("match_tcrs_to_db: loading database from ", db_tcrs_tsvfile)
    db_tcrs_df <- utils::read.delim(db_tcrs_tsvfile, stringsAsFactors = FALSE,
                                    row.names = NULL)
    db_tcrs_df <- .normalize_tcr_db_columns(db_tcrs_df)

    if (is.null(background_tcrs_df)) {
        background_tcrs_df <- tcr_df
    }

    find_significant_tcrdist_matches(
        query_tcrs_df = tcr_df[, c("va", "ja", "cdr3a", "cdr3a_nucseq",
                                    "vb", "jb", "cdr3b", "cdr3b_nucseq"),
                                drop = FALSE],
        db_tcrs_df = db_tcrs_df,
        organism = organism,
        adjusted_pvalue_threshold = adjusted_pvalue_threshold,
        background_tcrs_df = background_tcrs_df,
        num_random_samples = num_random_samples
    )
}


# ---------------------------------------------------------------------------
# strict_single_chain_match_tcrs_to_db
# ---------------------------------------------------------------------------

#' Strict CDR3 sequence matching against a database
#'
#' Finds database entries where the CDR3alpha or CDR3beta amino acid
#' sequence exactly matches a query TCR.
#'
#' @param tcr_df A data.frame with columns \code{cdr3a} and \code{cdr3b}.
#' @param organism Character string (e.g. \code{"human"}, \code{"mouse"}).
#' @param db_tcrs_tsvfile Path to a TSV file. If NULL, uses the built-in
#'   organism-specific single-chain database.
#' @return A list with elements \code{alpha_matches} and
#'   \code{beta_matches}, each a data.frame of matching database rows.
#' @export
strict_single_chain_match_tcrs_to_db <- function(tcr_df, organism,
                                                   db_tcrs_tsvfile = NULL) {
    if (is.null(db_tcrs_tsvfile)) {
        db_file <- if (organism == "human") {
            "human_tcr_db_for_matching.tsv"
        } else if (organism == "mouse") {
            "mouse_tcr_db_for_matching.tsv"
        } else {
            warning("strict_single_chain_match_tcrs_to_db: no default ",
                    "database for organism '", organism, "'")
            return(list(alpha_matches = data.frame(),
                        beta_matches = data.frame()))
        }
        db_tcrs_tsvfile <- system.file("extdata", db_file,
                                        package = "tcrdistR")
        if (!nzchar(db_tcrs_tsvfile)) {
            stop("strict_single_chain_match_tcrs_to_db: database file '",
                 db_file, "' not found. Ensure the package is properly ",
                 "installed.")
        }
        message("strict_single_chain_match_tcrs_to_db: matching to default ",
                organism, " database")
    }

    message("strict_single_chain_match_tcrs_to_db: loading database from ",
            db_tcrs_tsvfile)
    db_tcrs_df <- utils::read.delim(db_tcrs_tsvfile, stringsAsFactors = FALSE,
                                    row.names = NULL)
    db_tcrs_df <- .normalize_tcr_db_columns(db_tcrs_df)

    alpha_matches <- db_tcrs_df[db_tcrs_df$cdr3a %in% tcr_df$cdr3a, ,
                                 drop = FALSE]
    beta_matches <- db_tcrs_df[db_tcrs_df$cdr3b %in% tcr_df$cdr3b, ,
                                drop = FALSE]

    list(alpha_matches = alpha_matches, beta_matches = beta_matches)
}
