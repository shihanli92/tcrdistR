# TCR junction parsing, allele optimization, and background chain resampling.
#
# Extracted from rconga/R/tcr_sampler.R -- only functions needed for the
# clumping pipeline.  All computation-heavy paths delegate to C++ via
# rcpp_count_nuc_matches(), rcpp_find_alternate_alleles_batch(), and
# rcpp_resample_shuffled_tcr_chains().

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

.DEFAULT_MISMATCH_SCORE_JUNCTION_ANALYSIS <- -4L


# ---------------------------------------------------------------------------
# D-gene nucleotide sequence cache (built lazily)
# ---------------------------------------------------------------------------

#' Get TRBD nucleotide sequences for an organism
#'
#' Builds and caches a named list mapping 1-based integer D-gene IDs to their
#' nucleotide sequences.
#'
#' @param organism Character string. Organism identifier.
#' @return A named list with integer keys (as character names) and lowercase
#'   nucleotide sequence values.
#' @keywords internal
.get_trbd_nucseq <- function(organism) {
    cache_key <- paste0("trbd_nucseq_", organism)
    if (!is.null(.tcrdistR_env[[cache_key]])) {
        return(.tcrdistR_env[[cache_key]])
    }

    db <- load_gene_database(organism)
    d_ids <- sort(names(db)[vapply(db, function(g) {
        g$region == "D" && g$chain == "B"
    }, logical(1L))])

    result <- list()
    for (i in seq_along(d_ids)) {
        result[[as.character(i)]] <- db[[d_ids[i]]]$nucseq
    }

    .tcrdistR_env[[cache_key]] <- result
    result
}


# ---------------------------------------------------------------------------
# 1. .get_v_cdr3_nucseq
# ---------------------------------------------------------------------------

#' Get V-gene CDR3 region nucleotide sequence
#'
#' Returns the V-gene CDR3 region nucleotide sequence (lowercase), starting
#' from the conserved C codon.
#'
#' @param organism Character string. Organism identifier.
#' @param v_gene Character string. V-gene identifier including allele.
#' @return Character string. Lowercase nucleotide sequence starting at the
#'   conserved C codon.  Returns \code{""} if the gene is not found.
#' @keywords internal
.get_v_cdr3_nucseq <- function(organism, v_gene) {
    db <- load_gene_database(organism)
    vg <- db[[v_gene]]
    if (is.null(vg)) return("")

    v_nucseq        <- vg$nucseq
    v_nucseq_offset <- vg$nucseq_offset

    # Trim nucseq by offset (0-based offset -> 1-based substring start)
    v_nucseq <- substring(v_nucseq, v_nucseq_offset + 1L)

    v_alseq <- vg$alseq

    # cdr_columns is a list; last element's first value is the 1-based column
    # for the conserved C position in the alignment
    alseq_cpos <- vg$cdr_columns[[length(vg$cdr_columns)]][1L] - 1L  # 0-based

    # Count gaps before the C position in the alignment
    alseq_prefix <- substring(v_alseq, 1L, alseq_cpos)
    numgaps <- nchar(gsub("[^.]", "", alseq_prefix))
    v_cpos <- alseq_cpos - numgaps

    # Trim to start at C codon
    substring(v_nucseq, 3L * v_cpos + 1L)
}


# ---------------------------------------------------------------------------
# 2. .get_j_cdr3_nucseq
# ---------------------------------------------------------------------------

#' Get J-gene CDR3 region nucleotide sequence
#'
#' Returns the J-gene CDR3 region nucleotide sequence up to (but not
#' including) the GXG motif.
#'
#' @param organism Character string. Organism identifier.
#' @param j_gene Character string. J-gene identifier including allele.
#' @return Character string. Lowercase nucleotide sequence. Returns \code{""}
#'   if the gene is not found.
#' @keywords internal
.get_j_cdr3_nucseq <- function(organism, j_gene) {
    db <- load_gene_database(organism)
    jg <- db[[j_gene]]
    if (is.null(jg)) return("")

    j_nucseq        <- jg$nucseq
    j_nucseq_offset <- jg$nucseq_offset

    # Number of genome J amino acids in the CDR3 loop (up to but not
    # including GXG): count non-gap characters in the first CDR
    num_genome_j_aas_in_loop <- nchar(gsub("\\.", "", jg$cdrs[[1L]]))

    # Trim j_nucseq so that it extends up to the F/W position
    substring(j_nucseq, 1L, 3L * num_genome_j_aas_in_loop + j_nucseq_offset)
}


# ---------------------------------------------------------------------------
# 3. .analyze_junction
# ---------------------------------------------------------------------------

#' Analyze a TCR junction region
#'
#' Aligns the CDR3 nucleotide sequence against germline V and J sequences
#' to identify trimming, N-insertion lengths, and D-gene usage (beta only).
#'
#' @param organism Character string. Organism identifier.
#' @param v_gene Character string. V-gene identifier including allele.
#' @param j_gene Character string. J-gene identifier including allele.
#' @param cdr3_protseq Character string. CDR3 amino acid sequence.
#' @param cdr3_nucseq Character string. CDR3 nucleotide sequence (lowercase).
#' @param force_d_id Integer. If non-zero, only consider this D-gene ID.
#' @param mismatch_score Integer. Score for a mismatch.
#' @return A named list with elements: \code{new_nucseq},
#'   \code{cdr3_protseq_masked}, \code{cdr3_protseq_new_nuc_countstring},
#'   \code{cdr3_nucseq_src}, \code{trims}, \code{inserts}.
#' @keywords internal
.analyze_junction <- function(organism, v_gene, j_gene,
                               cdr3_protseq, cdr3_nucseq,
                               force_d_id = 0L,
                               mismatch_score = .DEFAULT_MISMATCH_SCORE_JUNCTION_ANALYSIS) {

    db <- load_gene_database(organism)
    ab <- db[[v_gene]]$chain

    v_nucseq <- .get_v_cdr3_nucseq(organism, v_gene)
    j_nucseq <- .get_j_cdr3_nucseq(organism, j_gene)

    cdr3_len <- nchar(cdr3_nucseq)

    # How far out do we match from V end
    num_matched_v <- rcpp_count_nuc_matches(v_nucseq, cdr3_nucseq,
                                            mismatch_score)

    # How far out do we match from J end (reverse alignment)
    j_rev <- paste0(rev(strsplit(j_nucseq, "", fixed = TRUE)[[1L]]),
                    collapse = "")
    cdr3_rev <- paste0(rev(strsplit(cdr3_nucseq, "", fixed = TRUE)[[1L]]),
                       collapse = "")
    num_matched_j <- rcpp_count_nuc_matches(j_rev, cdr3_rev, mismatch_score)

    # Handle overlap
    if (num_matched_v + num_matched_j > cdr3_len) {
        extra <- num_matched_v + num_matched_j - cdr3_len
        fake_v_trim <- extra %/% 2L
        fake_j_trim <- extra - fake_v_trim
        num_matched_v <- num_matched_v - fake_v_trim
        num_matched_j <- num_matched_j - fake_j_trim
    }

    stopifnot(num_matched_v + num_matched_j <= cdr3_len)

    if (num_matched_v + num_matched_j == cdr3_len) {
        nseq <- ""
    } else {
        nseq <- substring(cdr3_nucseq,
                          num_matched_v + 1L,
                          cdr3_len - num_matched_j)
    }

    # Build source annotation and new-nucleotide count vectors
    ncount <- rep(1L, cdr3_len)
    cdr3_nucseq_src <- rep("N", cdr3_len)

    if (num_matched_v > 0L) {
        for (i in seq_len(num_matched_v)) {
            ncount[i] <- 0L
            cdr3_nucseq_src[i] <- "V"
        }
    }
    if (num_matched_j > 0L) {
        for (i in seq_len(num_matched_j)) {
            ncount[cdr3_len + 1L - i] <- 0L
            cdr3_nucseq_src[cdr3_len + 1L - i] <- "J"
        }
    }

    v_trim <- nchar(v_nucseq) - num_matched_v
    j_trim <- nchar(j_nucseq) - num_matched_j

    # Initialize D/insertion variables
    n_vj_insert <- 0L
    n_vd_insert <- 0L
    n_dj_insert <- 0L
    d0_trim     <- 0L
    d1_trim     <- 0L
    best_d_id   <- 0L

    if (ab == "A") {
        n_vj_insert <- nchar(nseq)

    } else if (ab == "B") {
        # Look for D-gene segments in the N-region
        all_trbd <- .get_trbd_nucseq(organism)
        max_overlap <- 0L
        best_overlap_seq <- ""
        best_trim <- c(0L, 0L)
        nseq_len <- nchar(nseq)

        # Check if cdr3_nucseq has non-standard bases
        cdr3_nucseq_has_other <- grepl("[^acgt]", cdr3_nucseq)

        for (d_id_str in names(all_trbd)) {
            d_id <- as.integer(d_id_str)
            if (force_d_id != 0L && d_id != force_d_id) next

            d_nucseq <- all_trbd[[d_id_str]]
            d_len <- nchar(d_nucseq)

            for (start in seq_len(d_len) - 1L) {  # 0-based start
                for (stop in seq.int(start, d_len - 1L)) {  # 0-based stop
                    overlap_seq <- substring(d_nucseq, start + 1L, stop + 1L)
                    overlap_len <- stop - start + 1L

                    # Check if overlap_seq is found in nseq
                    if (!cdr3_nucseq_has_other && nseq_len > 0L &&
                        grepl(overlap_seq, nseq, fixed = TRUE)) {
                        if (overlap_len > max_overlap) {
                            max_overlap      <- overlap_len
                            best_d_id        <- d_id
                            best_overlap_seq <- overlap_seq
                            best_trim        <- c(start, d_len - 1L - stop)
                        }
                    }
                }
            }
        }

        if (max_overlap > 0L) {
            # Find position of D overlap in nseq
            pos <- regexpr(best_overlap_seq, nseq, fixed = TRUE)
            pos <- as.integer(pos) - 1L  # 0-based position in nseq

            # Mark D positions in the source annotation and ncount
            nseq_chars <- strsplit(nseq, "", fixed = TRUE)[[1L]]
            for (i in seq_len(max_overlap)) {
                idx_in_cdr3 <- pos + num_matched_v + i  # 1-based index
                ncount[idx_in_cdr3] <- 0L
                cdr3_nucseq_src[idx_in_cdr3] <- "D"
                # Replace with '+' in nseq
                nseq_chars[pos + i] <- "+"
            }
            nseq <- paste0(nseq_chars, collapse = "")

            n_vd_insert <- pos
            n_dj_insert <- nseq_len - pos - nchar(best_overlap_seq)
            d0_trim <- best_trim[1L]
            d1_trim <- best_trim[2L]
        } else {
            best_d_id   <- 0L
            n_vd_insert <- 0L
            n_dj_insert <- 0L
            n_vj_insert <- nchar(nseq)
            d0_trim     <- 0L
            d1_trim     <- 0L
        }
    }

    # Build masked protein sequence
    cdr3_protseq_masked <- ""
    cdr3_protseq_new_nuc_countstring <- ""

    if (nchar(cdr3_protseq) > 0L) {
        prot_chars <- strsplit(cdr3_protseq, "", fixed = TRUE)[[1L]]
        for (i in seq_along(prot_chars)) {
            nc <- sum(ncount[((i - 1L) * 3L + 1L):(i * 3L)])
            cdr3_protseq_new_nuc_countstring <- paste0(
                cdr3_protseq_new_nuc_countstring, nc)
            if (nc > 1L) {
                cdr3_protseq_masked <- paste0(cdr3_protseq_masked,
                                              prot_chars[i])
            } else if (nc == 1L) {
                cdr3_protseq_masked <- paste0(cdr3_protseq_masked,
                                              tolower(prot_chars[i]))
            } else {
                cdr3_protseq_masked <- paste0(cdr3_protseq_masked, "-")
            }
        }
    }

    list(
        new_nucseq                       = nseq,
        cdr3_protseq_masked              = cdr3_protseq_masked,
        cdr3_protseq_new_nuc_countstring = cdr3_protseq_new_nuc_countstring,
        cdr3_nucseq_src                  = paste0(cdr3_nucseq_src, collapse = ""),
        trims                            = c(v_trim, d0_trim, d1_trim, j_trim),
        inserts                          = c(best_d_id, n_vd_insert,
                                             n_dj_insert, n_vj_insert)
    )
}


# ---------------------------------------------------------------------------
# 4. .parse_tcr_junctions
# ---------------------------------------------------------------------------

#' Parse junction regions for a batch of TCRs
#'
#' Batch wrapper around \code{.analyze_junction}. Processes each TCR and
#' returns a data.frame with junction analysis results.
#'
#' @param organism Character string. Organism identifier.
#' @param tcr_df A data.frame with columns \code{va, ja, cdr3a, cdr3a_nucseq,
#'   vb, jb, cdr3b, cdr3b_nucseq}.
#' @return A data.frame with junction info columns needed for resampling.
#' @keywords internal
.parse_tcr_junctions <- function(organism, tcr_df) {
    n <- nrow(tcr_df)

    # Pre-allocate result vectors
    out_clone_index          <- integer(n)
    out_va                   <- character(n)
    out_ja                   <- character(n)
    out_cdr3a                <- character(n)
    out_cdr3a_nucseq         <- character(n)
    out_cdr3a_protseq_masked <- character(n)
    out_cdr3a_nucseq_src     <- character(n)
    out_va_trim              <- integer(n)
    out_ja_trim              <- integer(n)
    out_a_insert             <- vector("list", n)
    out_vb                   <- character(n)
    out_jb                   <- character(n)
    out_cdr3b                <- character(n)
    out_cdr3b_nucseq         <- character(n)
    out_cdr3b_protseq_masked <- character(n)
    out_cdr3b_nucseq_src     <- character(n)
    out_vb_trim              <- integer(n)
    out_d0_trim              <- integer(n)
    out_d1_trim              <- integer(n)
    out_jb_trim              <- integer(n)
    out_vd_insert            <- vector("list", n)
    out_dj_insert            <- vector("list", n)
    out_vj_insert            <- vector("list", n)

    for (ii in seq_len(n)) {
        if ((ii - 1L) %% 1000L == 0L) {
            message(".parse_tcr_junctions: ", ii - 1L, " ", n)
        }

        va           <- tcr_df$va[ii]
        ja           <- tcr_df$ja[ii]
        cdr3a        <- tcr_df$cdr3a[ii]
        cdr3a_nucseq <- tcr_df$cdr3a_nucseq[ii]

        vb           <- tcr_df$vb[ii]
        jb           <- tcr_df$jb[ii]
        cdr3b        <- tcr_df$cdr3b[ii]
        cdr3b_nucseq <- tcr_df$cdr3b_nucseq[ii]

        aresults <- .analyze_junction(organism, va, ja, cdr3a, cdr3a_nucseq)
        bresults <- .analyze_junction(organism, vb, jb, cdr3b, cdr3b_nucseq)

        out_clone_index[ii]          <- ii - 1L
        out_va[ii]                   <- va
        out_ja[ii]                   <- ja
        out_cdr3a[ii]                <- cdr3a
        out_cdr3a_nucseq[ii]         <- cdr3a_nucseq
        out_cdr3a_protseq_masked[ii] <- aresults$cdr3_protseq_masked
        out_cdr3a_nucseq_src[ii]     <- aresults$cdr3_nucseq_src
        out_va_trim[ii]              <- aresults$trims[1L]
        out_ja_trim[ii]              <- aresults$trims[4L]
        out_a_insert[[ii]]           <- aresults$inserts[4L]
        out_vb[ii]                   <- vb
        out_jb[ii]                   <- jb
        out_cdr3b[ii]                <- cdr3b
        out_cdr3b_nucseq[ii]         <- cdr3b_nucseq
        out_cdr3b_protseq_masked[ii] <- bresults$cdr3_protseq_masked
        out_cdr3b_nucseq_src[ii]     <- bresults$cdr3_nucseq_src
        out_vb_trim[ii]              <- bresults$trims[1L]
        out_d0_trim[ii]              <- bresults$trims[2L]
        out_d1_trim[ii]              <- bresults$trims[3L]
        out_jb_trim[ii]              <- bresults$trims[4L]
        out_vd_insert[[ii]]          <- bresults$inserts[2L]
        out_dj_insert[[ii]]          <- bresults$inserts[3L]
        out_vj_insert[[ii]]          <- bresults$inserts[4L]
    }

    data.frame(
        clone_index          = out_clone_index,
        va                   = out_va,
        ja                   = out_ja,
        cdr3a                = out_cdr3a,
        cdr3a_nucseq         = out_cdr3a_nucseq,
        cdr3a_protseq_masked = out_cdr3a_protseq_masked,
        cdr3a_nucseq_src     = out_cdr3a_nucseq_src,
        va_trim              = out_va_trim,
        ja_trim              = out_ja_trim,
        a_insert             = unlist(out_a_insert),
        vb                   = out_vb,
        jb                   = out_jb,
        cdr3b                = out_cdr3b,
        cdr3b_nucseq         = out_cdr3b_nucseq,
        cdr3b_protseq_masked = out_cdr3b_protseq_masked,
        cdr3b_nucseq_src     = out_cdr3b_nucseq_src,
        vb_trim              = out_vb_trim,
        d0_trim              = out_d0_trim,
        d1_trim              = out_d1_trim,
        jb_trim              = out_jb_trim,
        vd_insert            = unlist(out_vd_insert),
        dj_insert            = unlist(out_dj_insert),
        vj_insert            = unlist(out_vj_insert),
        stringsAsFactors     = FALSE
    )
}


# ---------------------------------------------------------------------------
# 5. .find_alternate_alleles_for_tcrs
# ---------------------------------------------------------------------------

#' Find alternate alleles for a batch of TCRs
#'
#' R wrapper for \code{rcpp_find_alternate_alleles_batch}. For each TCR,
#' tries alternate alleles. If an alternate allele is consistently better
#' (count >= \code{min_better_count} AND ratio >= \code{min_better_ratio}:1),
#' swaps it into the data.
#'
#' @param organism Character string. Organism identifier.
#' @param tcr_df A data.frame with columns \code{va, ja, cdr3a_nucseq,
#'   vb, jb, cdr3b_nucseq}.
#' @param min_better_ratio Numeric. Minimum ratio for the alternate allele
#'   to be allowed.
#' @param min_better_count Integer. Minimum count for the alternate allele
#'   to be allowed.
#' @param min_improvement Integer. Minimum improvement in match count.
#' @param verbose Logical. If \code{TRUE}, print progress information.
#' @return A data.frame of the same structure as \code{tcr_df} with
#'   potentially updated gene names.
#' @keywords internal
.find_alternate_alleles_for_tcrs <- function(organism, tcr_df,
                                              min_better_ratio = 10,
                                              min_better_count = 5L,
                                              min_improvement = 2L,
                                              verbose = FALSE) {
    if (verbose) {
        message(".find_alternate_alleles_for_tcrs: num_tcrs: ", nrow(tcr_df))
    }

    n <- nrow(tcr_df)
    db <- load_gene_database(organism)
    all_names <- names(db)

    # Pre-compute V and J CDR3 nucseqs for all genes in database
    v_nucseqs <- vapply(all_names, function(g)
        tryCatch(.get_v_cdr3_nucseq(organism, g),
                 error = function(e) ""), character(1L))
    j_nucseqs <- vapply(all_names, function(g)
        tryCatch(.get_j_cdr3_nucseq(organism, g),
                 error = function(e) ""), character(1L))

    mismatch <- .DEFAULT_MISMATCH_SCORE_JUNCTION_ANALYSIS
    rcpp_result <- rcpp_find_alternate_alleles_batch(
        tcr_df$va, tcr_df$ja, tcr_df$cdr3a_nucseq,
        tcr_df$vb, tcr_df$jb, tcr_df$cdr3b_nucseq,
        all_names, v_nucseqs, j_nucseqs,
        mismatch, min_improvement
    )

    all_new_va <- as.character(rcpp_result$new_va)
    all_new_ja <- as.character(rcpp_result$new_ja)
    all_new_vb <- as.character(rcpp_result$new_vb)
    all_new_jb <- as.character(rcpp_result$new_jb)

    # Reconstruct all_counts nested list from C++ counts data.frame
    counts_df <- rcpp_result$counts
    all_counts <- list()
    for (row_i in seq_len(nrow(counts_df))) {
        old_g <- counts_df$old_gene[row_i]
        new_g <- counts_df$new_gene[row_i]
        cnt   <- counts_df$count[row_i]
        if (is.null(all_counts[[old_g]])) {
            all_counts[[old_g]] <- list()
        }
        all_counts[[old_g]][[new_g]] <- cnt
    }

    # Determine allowed swaps
    allowed_swaps <- character(0L)
    for (g in names(all_counts)) {
        counts_g <- all_counts[[g]]
        count_vec <- unlist(counts_g)
        sorted_idx <- order(count_vec, decreasing = TRUE)
        sorted_names <- names(count_vec)[sorted_idx]
        sorted_vals  <- count_vec[sorted_idx]

        if (sorted_names[1L] != g ||
            (length(sorted_names) > 1L &&
             sorted_vals[2L] >= sorted_vals[1L] %/% min_better_ratio &&
             sorted_vals[2L] >= min_better_count)) {
            if (sorted_names[1L] != g) {
                alt_g <- sorted_names[1L]
            } else {
                alt_g <- sorted_names[2L]
            }
            if (verbose) {
                message(".find_alternate_alleles_for_tcrs: ",
                        "allow alternate gene: ", g, " ", alt_g)
            }
            allowed_swaps <- c(allowed_swaps, alt_g)
        }
    }

    # Apply allowed swaps
    result <- tcr_df
    for (ii in seq_len(n)) {
        if (all_new_va[ii] %in% allowed_swaps) result$va[ii] <- all_new_va[ii]
        if (all_new_ja[ii] %in% allowed_swaps) result$ja[ii] <- all_new_ja[ii]
        if (all_new_vb[ii] %in% allowed_swaps) result$vb[ii] <- all_new_vb[ii]
        if (all_new_jb[ii] %in% allowed_swaps) result$jb[ii] <- all_new_jb[ii]
    }

    result
}


# ---------------------------------------------------------------------------
# 6. .resample_shuffled_tcr_chains
# ---------------------------------------------------------------------------

#' Resample shuffled TCR chains
#'
#' R wrapper for \code{rcpp_resample_shuffled_tcr_chains}. Generates random
#' TCR chains by recombining junction segments from existing TCRs.
#'
#' @param organism Character string. Organism identifier.
#' @param num_samples Integer. Number of resampled chains to generate.
#' @param chain Character string. Either \code{"A"} or \code{"B"}.
#' @param junctions_df A data.frame as produced by
#'   \code{.parse_tcr_junctions}.
#' @param preserve_vj_pairings Logical. If \code{TRUE}, sample chimeras
#'   preserving V-gene and J-gene pairings.
#' @return A data.frame with columns \code{v_gene, j_gene, cdr3,
#'   cdr3_nucseq}.
#' @keywords internal
.resample_shuffled_tcr_chains <- function(organism, num_samples, chain,
                                           junctions_df,
                                           preserve_vj_pairings = FALSE,
                                           max_attempts = 100L * num_samples) {
    stopifnot(chain %in% c("A", "B"))

    # Build junctions list: each element has v_gene, j_gene, nucseq,
    # breakpoints_pre_d, breakpoints_post_d
    junctions <- vector("list", nrow(junctions_df))

    for (row_idx in seq_len(nrow(junctions_df))) {
        l <- junctions_df[row_idx, , drop = FALSE]

        if (chain == "A") {
            nucseq_src <- l$cdr3a_nucseq_src
            vg <- l$va
            jg <- l$ja
            nucseq <- l$cdr3a_nucseq
        } else {
            nucseq_src <- l$cdr3b_nucseq_src
            vg <- l$vb
            jg <- l$jb
            nucseq <- l$cdr3b_nucseq
        }

        src_chars <- strsplit(nucseq_src, "", fixed = TRUE)[[1L]]
        src_len <- length(src_chars)

        breakpoints_pre_d  <- integer(0L)
        breakpoints_post_d <- integer(0L)
        post_d <- FALSE

        for (ii in seq.int(2L, src_len)) {
            a <- src_chars[ii - 1L]
            b <- src_chars[ii]
            if (a == "D") post_d <- TRUE
            if (a == b && a != "N") {
                next
            }
            # Store Python-style 0-based breakpoint indices
            ii_py <- ii - 1L
            bp_pos <- ii_py
            bp_neg <- ii_py - src_len

            if (post_d) {
                breakpoints_post_d <- c(breakpoints_post_d, bp_pos, bp_neg)
            } else {
                breakpoints_pre_d <- c(breakpoints_pre_d, bp_pos, bp_neg)
            }
        }

        # Remove duplicates
        breakpoints_pre_d  <- unique(breakpoints_pre_d)
        breakpoints_post_d <- unique(breakpoints_post_d)

        junctions[[row_idx]] <- list(
            v_gene             = vg,
            j_gene             = jg,
            nucseq             = nucseq,
            breakpoints_pre_d  = breakpoints_pre_d,
            breakpoints_post_d = breakpoints_post_d
        )
    }

    # ---- C++ fast path (non-preserve_vj_pairings only) ----
    if (!preserve_vj_pairings) {
        j_v_genes <- vapply(junctions, `[[`, character(1L), "v_gene")
        j_j_genes <- vapply(junctions, `[[`, character(1L), "j_gene")
        j_nucseqs <- vapply(junctions, `[[`, character(1L), "nucseq")
        j_bp_pre  <- lapply(junctions, `[[`, "breakpoints_pre_d")
        j_bp_post <- lapply(junctions, `[[`, "breakpoints_post_d")

        return(rcpp_resample_shuffled_tcr_chains(
            j_v_genes, j_j_genes, j_nucseqs,
            j_bp_pre, j_bp_post,
            chain, num_samples,
            max_attempts
        ))
    }

    # ---- Pure-R fallback for preserve_vj_pairings = TRUE ----
    message(".resample_shuffled_tcr_chains: create gene --> junctions mapping")
    v_gene2junctions <- list()
    j_gene2junctions <- list()

    for (j in junctions) {
        vg_trimmed <- sub("\\*.*$", "", j$v_gene)
        jg_trimmed <- sub("\\*.*$", "", j$j_gene)

        if (is.null(v_gene2junctions[[vg_trimmed]])) {
            v_gene2junctions[[vg_trimmed]] <- list()
        }
        v_gene2junctions[[vg_trimmed]] <- c(
            v_gene2junctions[[vg_trimmed]], list(j))

        if (is.null(j_gene2junctions[[jg_trimmed]])) {
            j_gene2junctions[[jg_trimmed]] <- list()
        }
        j_gene2junctions[[jg_trimmed]] <- c(
            j_gene2junctions[[jg_trimmed]], list(j))
    }

    n_junctions <- length(junctions)
    new_tcrs <- vector("list", num_samples)
    count <- 0L
    attempts <- 0L
    successes <- 0L

    while (count < num_samples && attempts < max_attempts) {
        j <- junctions[[sample.int(n_junctions, 1L)]]
        vg_trimmed <- sub("\\*.*$", "", j$v_gene)
        jg_trimmed <- sub("\\*.*$", "", j$j_gene)
        v_pool <- v_gene2junctions[[vg_trimmed]]
        j_pool <- j_gene2junctions[[jg_trimmed]]
        t1 <- v_pool[[sample.int(length(v_pool), 1L)]]
        t2 <- j_pool[[sample.int(length(j_pool), 1L)]]

        nucseq1 <- t1$nucseq
        nucseq2 <- t2$nucseq
        if (nucseq1 == nucseq2) { attempts <- attempts + 1L; next }

        # Choose breakpoint indices to try
        if (chain == "A") {
            inds <- 1L  # only pre_d
        } else {
            inds <- sample(c(1L, 2L))  # randomize pre_d vs post_d
        }

        success <- FALSE
        for (ind in inds) {
            if (ind == 1L) {
                bp_set1 <- t1$breakpoints_pre_d
                bp_set2 <- t2$breakpoints_pre_d
            } else {
                bp_set1 <- t1$breakpoints_post_d
                bp_set2 <- t2$breakpoints_post_d
            }

            shared <- intersect(bp_set1, bp_set2)
            if (length(shared) > 0L) {
                bp <- shared[sample.int(length(shared), 1L)]

                # Construct chimeric nucleotide sequence (Python-style bp)
                if (bp >= 0L) {
                    nucseq <- paste0(
                        substring(nucseq1, 1L, bp),
                        substring(nucseq2, bp + 1L))
                } else {
                    len1 <- nchar(nucseq1)
                    len2 <- nchar(nucseq2)
                    cut1 <- len1 + bp
                    cut2 <- len2 + bp
                    nucseq <- paste0(
                        substring(nucseq1, 1L, cut1),
                        substring(nucseq2, cut2 + 1L))
                }

                # Check for multiples of 3
                if (nchar(nucseq) %% 3L != 0L) next

                # Translate and check for stop codons
                cdr3 <- get_translation(nucseq)
                if (!grepl("\\*", cdr3)) {
                    success <- TRUE
                    count <- count + 1L
                    new_tcrs[[count]] <- data.frame(
                        v_gene      = t1$v_gene,
                        j_gene      = t2$j_gene,
                        cdr3        = cdr3,
                        cdr3_nucseq = nucseq,
                        stringsAsFactors = FALSE
                    )
                    break
                }
            }
        }

        attempts <- attempts + 1L
        successes <- successes + as.integer(success)
        if (success && count >= num_samples) break
    }

    if (count < num_samples) {
        warning(sprintf(
            ".resample_shuffled_tcr_chains: reached max_attempts (%d) with only %d/%d samples (success_rate=%.2f%%)",
            max_attempts, count, num_samples,
            if (attempts > 0L) 100.0 * successes / attempts else 0.0),
            call. = FALSE)
    }

    message(sprintf(".resample_shuffled_tcr_chains: success_rate: %.2f",
                    if (attempts > 0L) 100.0 * successes / attempts else 0.0))
    do.call(rbind, new_tcrs[seq_len(count)])
}
