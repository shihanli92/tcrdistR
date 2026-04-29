# Ported from Python CoNGA: conga/tcrdist/all_genes.py

# Internal separator used in the TSV cdr_columns / cdrs fields.
.CDRS_SEP <- ";"

# ---------------------------------------------------------------------------
# Helper: parse one TSV row into a gene entry list
# ---------------------------------------------------------------------------

#' @keywords internal
.parse_gene_row <- function(id, organism, chain, region, nucseq, frame,
                            aligned_protseq, cdr_columns_str, cdrs_str) {
    # ---- CDRs and CDR column ranges ----------------------------------------
    if (is.na(cdrs_str) || cdrs_str == "") {
        cdrs        <- character(0L)
        cdr_columns <- list()
    } else {
        cdrs <- strsplit(cdrs_str, .CDRS_SEP, fixed = TRUE)[[1L]]
        cdr_columns <- lapply(
            strsplit(cdr_columns_str, .CDRS_SEP, fixed = TRUE)[[1L]],
            function(x) as.integer(strsplit(x, "-", fixed = TRUE)[[1L]])
        )
    }

    # ---- nucleotide offset (0-based, matching Python nucseq_offset) --------
    nucseq_offset <- as.integer(frame) - 1L

    # ---- protein translation -----------------------------------------------
    frame_str <- paste0("+", as.integer(frame))
    protseq   <- get_translation(nucseq, frame_str)

    # ---- assertions (abort loudly to catch data / porting issues) ----------
    expected_protseq <- gsub(".", "", aligned_protseq, fixed = TRUE)
    if (protseq != expected_protseq) {
        stop(sprintf(
            "Translation mismatch for gene '%s': got '%s', expected '%s'",
            id, protseq, expected_protseq
        ))
    }
    if (length(cdrs) > 0L) {
        for (i in seq_along(cdrs)) {
            extracted <- substr(aligned_protseq,
                                cdr_columns[[i]][1L],
                                cdr_columns[[i]][2L])
            if (extracted != cdrs[[i]]) {
                stop(sprintf(
                    "CDR mismatch for gene '%s' CDR%d: got '%s', expected '%s'",
                    id, i, extracted, cdrs[[i]]
                ))
            }
        }
    }

    list(
        id            = id,
        organism      = organism,
        chain         = chain,
        region        = region,
        nucseq        = nucseq,
        alseq         = aligned_protseq,
        cdrs          = cdrs,
        cdr_columns   = cdr_columns,
        nucseq_offset = nucseq_offset,
        protseq       = protseq,
        rep           = NA_character_,
        mm1_rep       = NA_character_,
        count_rep     = NA_character_
    )
}

# ---------------------------------------------------------------------------
# Helper: compute representatives for V genes in one organism+chain
# ---------------------------------------------------------------------------

#' @keywords internal
.compute_v_reps <- function(genes_vc) {
    # genes_vc: named list of gene-entry lists, all region == 'V', same chain

    ids <- names(genes_vc)
    if (length(ids) == 0L) return(genes_vc)

    # ---- Build merged loop sequences (all CDRs except the last, space-joined)
    merged_loopseqs <- vapply(ids, function(id) {
        cdrs <- genes_vc[[id]]$cdrs
        if (length(cdrs) <= 1L) return("")
        paste(cdrs[-length(cdrs)], collapse = " ")
    }, character(1L))

    # ---- Exact neighbours and mm1 neighbours --------------------------------
    all_loopseq_nbrs     <- vector("list", length(ids))
    all_loopseq_nbrs_mm1 <- vector("list", length(ids))
    names(all_loopseq_nbrs)     <- ids
    names(all_loopseq_nbrs_mm1) <- ids

    for (id1 in ids) {
        seq1   <- merged_loopseqs[[id1]]
        g1     <- genes_vc[[id1]]
        # cpos: 0-indexed start of last CDR (Python: cdr_columns[-1][0] - 1)
        # used only to bound the all_mismatches slice: alseq[:cpos+2] in Python
        # = substr(alseq, 1, cpos + 2) in R  where cpos = last_cdr_start_R - 1
        cpos   <- g1$cdr_columns[[length(g1$cdr_columns)]][1L] - 1L
        alseq1 <- g1$alseq

        chars1 <- strsplit(seq1, "", fixed = TRUE)[[1L]]
        nchar1 <- length(chars1)

        exact_nbrs <- character(0)
        mm1_nbrs   <- character(0)

        for (id2 in ids) {
            seq2 <- merged_loopseqs[[id2]]

            if (seq1 == seq2) {
                exact_nbrs <- c(exact_nbrs, id2)
                mm1_nbrs   <- c(mm1_nbrs,   id2)
                next
            }

            # ---- Count loop-sequence mismatches ----------------------------
            chars2 <- strsplit(seq2, "", fixed = TRUE)[[1L]]
            nchar2 <- length(chars2)
            stopifnot(nchar1 == nchar2)  # must be equal (same organism/chain)
            compare_len <- nchar1

            loop_mismatches      <- 0L
            loop_mismatches_cdrx <- 0L
            loop_mismatch_seqs   <- list()
            spaces               <- 0L
            too_many             <- FALSE

            for (i in seq_len(compare_len)) {
                a <- chars1[i]
                b <- chars2[i]
                if (a == " ") {
                    spaces <- spaces + 1L
                    next
                }
                if (a != b) {
                    # Gaps or stop codons count heavily
                    if (grepl("[*.]", a, fixed = FALSE) ||
                        grepl("[*.]", b, fixed = FALSE)) {
                        loop_mismatches <- loop_mismatches + 10L
                        too_many <- TRUE
                        break
                    }
                    if (spaces <= 1L) {
                        # CDR1 or CDR2 mismatch
                        loop_mismatches <- loop_mismatches + 1L
                        loop_mismatch_seqs <- c(loop_mismatch_seqs, list(c(a, b)))
                    } else {
                        # CDR3 stub mismatch
                        loop_mismatches_cdrx <- loop_mismatches_cdrx + 1L
                    }
                    if (loop_mismatches > 1L) {
                        too_many <- TRUE
                        break
                    }
                }
            }

            if (too_many || loop_mismatches > 1L) next

            # ---- Count all mismatches up to CDR3+2 boundary ----------------
            alseq2        <- genes_vc[[id2]]$alseq
            slice_end     <- cpos + 2L  # = last_cdr_start_R + 1  (1-indexed inclusive)
            alseq1_slice  <- substr(alseq1, 1L, slice_end)
            alseq2_slice  <- substr(alseq2, 1L, slice_end)
            achars1       <- strsplit(alseq1_slice, "", fixed = TRUE)[[1L]]
            achars2       <- strsplit(alseq2_slice, "", fixed = TRUE)[[1L]]
            compare_alen  <- min(length(achars1), length(achars2))

            all_mismatches <- 0L
            for (i in seq_len(compare_alen)) {
                if (achars1[i] != achars2[i]) {
                    if (grepl("[*.]", achars1[i], fixed = FALSE) ||
                        grepl("[*.]", achars2[i], fixed = FALSE)) {
                        all_mismatches <- all_mismatches + 10L
                    } else {
                        all_mismatches <- all_mismatches + 1L
                    }
                }
            }

            # ---- Apply mm1 acceptance criteria -----------------------------
            if (loop_mismatches <= 1L &&
                loop_mismatches + loop_mismatches_cdrx <= 2L &&
                all_mismatches <= 10L) {

                blscore <- if (loop_mismatches == 1L) {
                    pair <- loop_mismatch_seqs[[1L]]
                    rcpp_blosum62_lookup(pair[1L], pair[2L])
                } else {
                    100L
                }

                if (blscore >= 1L) {
                    mm1_nbrs <- c(mm1_nbrs, id2)
                }
            }
        }  # end id2 loop

        all_loopseq_nbrs[[id1]]     <- exact_nbrs
        all_loopseq_nbrs_mm1[[id1]] <- mm1_nbrs

    }  # end id1 loop

    # ---- Transitive closure of mm1 neighbours -------------------------------
    # Mirrors the Python 'while True' loop with one-neighbour-at-a-time expansion
    repeat {
        new_nbrs <- FALSE
        for (id1 in ids) {
            new_id1_nbrs <- FALSE
            for (id2 in all_loopseq_nbrs_mm1[[id1]]) {
                for (id3 in all_loopseq_nbrs_mm1[[id2]]) {
                    if (!(id3 %in% all_loopseq_nbrs_mm1[[id1]])) {
                        all_loopseq_nbrs_mm1[[id1]] <-
                            c(all_loopseq_nbrs_mm1[[id1]], id3)
                        new_id1_nbrs <- TRUE
                        break
                    }
                }
                if (new_id1_nbrs) break
            }
            if (new_id1_nbrs) new_nbrs <- TRUE
        }
        if (!new_nbrs) break
    }

    # ---- Assign rep and mm1_rep (byte-order minimum of neighbour set) ------
    # Use sort(method = "radix") to match Python's byte-order min(), since
    # R's default min() uses locale collation which differs for punctuation.
    for (id in ids) {
        genes_vc[[id]]$rep     <- sort(all_loopseq_nbrs[[id]], method = "radix")[1L]
        genes_vc[[id]]$mm1_rep <- sort(all_loopseq_nbrs_mm1[[id]], method = "radix")[1L]
    }

    genes_vc
}

# ---------------------------------------------------------------------------
# Helper: compute representatives for J genes in one organism+chain
# ---------------------------------------------------------------------------

#' @keywords internal
.compute_j_reps <- function(genes_jc) {
    # genes_jc: named list of gene-entry lists, all region == 'J', same chain

    ids <- names(genes_jc)
    if (length(ids) == 0L) return(genes_jc)

    # ---- Build J loop sequences: protseq[1 .. num+3] -----------------------
    # num = number of non-gap residues in the sole CDR
    # The +3 extension captures the conserved GXG motif that follows the CDR
    jloopseqs <- vapply(ids, function(id) {
        g <- genes_jc[[id]]
        if (length(g$cdrs) == 0L) return("")
        num <- nchar(gsub(".", "", g$cdrs[1L], fixed = TRUE))
        substr(g$protseq, 1L, num + 3L)
    }, character(1L))

    # ---- Exact neighbours ---------------------------------------------------
    all_jloopseq_nbrs     <- vector("list", length(ids))
    all_jloopseq_nbrs_mm1 <- vector("list", length(ids))
    names(all_jloopseq_nbrs)     <- ids
    names(all_jloopseq_nbrs_mm1) <- ids

    for (id1 in ids) {
        seq1       <- jloopseqs[[id1]]
        exact_nbrs <- character(0)
        mm1_nbrs   <- character(0)

        for (id2 in ids) {
            seq2 <- jloopseqs[[id2]]
            if (seq1 == seq2) {
                exact_nbrs <- c(exact_nbrs, id2)
                mm1_nbrs   <- c(mm1_nbrs,   id2)
            }
        }
        all_jloopseq_nbrs[[id1]]     <- exact_nbrs
        all_jloopseq_nbrs_mm1[[id1]] <- mm1_nbrs
    }

    # ---- Assign rep and mm1_rep (byte-order minimum) -------------------------
    for (id in ids) {
        genes_jc[[id]]$rep     <- sort(all_jloopseq_nbrs[[id]], method = "radix")[1L]
        genes_jc[[id]]$mm1_rep <- sort(all_jloopseq_nbrs_mm1[[id]], method = "radix")[1L]
    }

    genes_jc
}

# ---------------------------------------------------------------------------
# Helper: build the full gene database (called once, result cached)
# ---------------------------------------------------------------------------

#' @keywords internal
.build_gene_database <- function() {
    db_path <- system.file(
        "extdata", "combo_xcr_2023-12-30.tsv",
        package = "tcrdistR"
    )
    # Fallback for development mode (devtools::load_all or direct sourcing)
    if (!nzchar(db_path)) {
        dev_path <- file.path(
            getOption("tcrdistR.dev_root",
                      default = getwd()),
            "inst", "extdata", "combo_xcr_2023-12-30.tsv"
        )
        if (file.exists(dev_path)) {
            db_path <- dev_path
        } else {
            stop("Gene database file 'combo_xcr_2023-12-30.tsv' not found. ",
                 "Install the package or set options(tcrdistR.dev_root = '/path/to/tcrdistR').")
        }
    }

    dt <- utils::read.delim(db_path, stringsAsFactors = FALSE)

    # ---- Parse all genes into nested list: all_genes[[organism]][[id]] -----
    all_genes <- list()

    for (i in seq_len(nrow(dt))) {
        row <- dt[i, , drop = FALSE]
        g   <- .parse_gene_row(
            id              = row$id,
            organism        = row$organism,
            chain           = row$chain,
            region          = row$region,
            nucseq          = row$nucseq,
            frame           = row$frame,
            aligned_protseq = row$aligned_protseq,
            cdr_columns_str = row$cdr_columns,
            cdrs_str        = row$cdrs
        )
        org <- g$organism
        if (is.null(all_genes[[org]])) {
            all_genes[[org]] <- list()
        }
        all_genes[[org]][[g$id]] <- g
    }

    # ---- Compute representatives per organism and chain --------------------
    organisms <- names(all_genes)

    for (org in organisms) {
        genes_org <- all_genes[[org]]
        ids_org   <- names(genes_org)

        for (ch in c("A", "B")) {
            # --- V genes ---
            v_ids <- ids_org[vapply(ids_org,
                                    function(id) genes_org[[id]]$chain == ch &&
                                                 genes_org[[id]]$region == "V",
                                    logical(1L))]
            if (length(v_ids) > 0L) {
                genes_vc <- genes_org[v_ids]
                genes_vc <- .compute_v_reps(genes_vc)
                for (id in v_ids) {
                    all_genes[[org]][[id]] <- genes_vc[[id]]
                }
            }

            # --- J genes ---
            j_ids <- ids_org[vapply(ids_org,
                                    function(id) genes_org[[id]]$chain == ch &&
                                                 genes_org[[id]]$region == "J",
                                    logical(1L))]
            if (length(j_ids) > 0L) {
                genes_jc <- genes_org[j_ids]
                genes_jc <- .compute_j_reps(genes_jc)
                for (id in j_ids) {
                    all_genes[[org]][[id]] <- genes_jc[[id]]
                }
            }

            # --- D genes: no loopseq neighbours; set rep = mm1_rep = id ----
            d_ids <- ids_org[vapply(ids_org,
                                    function(id) genes_org[[id]]$chain == ch &&
                                                 genes_org[[id]]$region == "D",
                                    logical(1L))]
            for (id in d_ids) {
                all_genes[[org]][[id]]$rep     <- id
                all_genes[[org]][[id]]$mm1_rep <- id
            }
        }

        # ---- count_rep: always trim_allele_to_gene(id) (CLASSIC_COUNTREPS=FALSE)
        for (id in ids_org) {
            all_genes[[org]][[id]]$count_rep <- trim_allele_to_gene(id)
        }
    }

    all_genes
}

# ---------------------------------------------------------------------------
# Internal: compute V-region distance matrix for one organism+chain
# ---------------------------------------------------------------------------

#' @keywords internal
.compute_v_region_distance_matrix <- function(organism, chain) {
    cache_key <- paste0("v_dist_", organism, "_", chain)

    # Return cached result if available
    cached <- .tcrdistR_env$v_dist_matrices[[cache_key]]
    if (!is.null(cached)) {
        return(cached)
    }

    # Load the gene database for this organism
    all_genes_org <- load_gene_database(organism)

    # Filter V genes for this chain
    ids <- names(all_genes_org)
    v_ids <- ids[vapply(ids, function(id) {
        g <- all_genes_org[[id]]
        g$chain == chain && g$region == "V"
    }, logical(1L))]

    if (length(v_ids) == 0L) {
        stop(sprintf(
            ".compute_v_region_distance_matrix: no V genes found for organism '%s', chain '%s'",
            organism, chain
        ))
    }

    # Build merged loop sequences: all CDRs except the last, joined by space
    # This matches compute_all_v_region_distances in rconga/R/tcr_distances.R
    loop_seqs <- vapply(v_ids, function(id) {
        cdrs <- all_genes_org[[id]]$cdrs
        if (length(cdrs) > 1L) {
            paste(cdrs[-length(cdrs)], collapse = " ")
        } else {
            ""
        }
    }, character(1L))

    # Call C++ function to compute the distance matrix
    result <- rcpp_compute_v_region_distances(
        v_ids,
        loop_seqs,
        as.numeric(WEIGHT_V_REGION),
        as.numeric(GAP_PENALTY_V_REGION)
    )

    # Cache and return
    .tcrdistR_env$v_dist_matrices[[cache_key]] <- result
    result
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Trim allele designation from a gene ID
#'
#' Removes the allele suffix (everything from \code{*} onward) from a TCR or
#' BCR gene identifier, returning only the gene-level name.
#'
#' @param gene_id Character string. A gene identifier containing an allele
#'   designation, e.g. \code{"TRAV1*01"} or \code{"IGHV1-18*03"}.
#' @return A character string with the allele suffix stripped, e.g.
#'   \code{"TRAV1"} or \code{"IGHV1-18"}.
#' @examples
#' trim_allele_to_gene("TRAV1*01")    # "TRAV1"
#' trim_allele_to_gene("IGHV1-18*03") # "IGHV1-18"
#' @seealso \code{\link{load_gene_database}}
#' @export
trim_allele_to_gene <- function(gene_id) {
    sub("\\*.*$", "", gene_id)
}

#' Load and cache the TCR/BCR gene database
#'
#' Parses the bundled gene database TSV on first call, computes sequence
#' representatives (exact and mm1) for V and J genes, and caches the result
#' in the package-private environment \code{.tcrdistR_env$all_genes}. Subsequent
#' calls return the cached data without re-parsing.
#'
#' The database file is \code{inst/extdata/combo_xcr_2023-12-30.tsv} and
#' contains 2836 gene entries across organisms including \code{"human"},
#' \code{"mouse"}, \code{"human_ig"}, \code{"mouse_ig"}, \code{"human_gd"},
#' \code{"mouse_gd"}, and \code{"rhesus"}.
#'
#' @param organism Character string or \code{NULL}. If a non-\code{NULL}
#'   string is supplied, only the gene entries for that organism are returned.
#'   If \code{NULL} (default), the complete nested list for all organisms is
#'   returned.
#' @return
#' When \code{organism} is \code{NULL}: a named list keyed by organism name,
#' where each element is itself a named list of gene entries keyed by gene ID.
#'
#' When \code{organism} is a character string: a named list of gene entries for
#' that organism, keyed by gene ID.
#'
#' Each gene entry is a named list with fields:
#' \describe{
#'   \item{id}{Character. Gene identifier including allele, e.g. \code{"TRAV1*01"}.}
#'   \item{organism}{Character. Organism name.}
#'   \item{chain}{Character. Either \code{"A"} (alpha / gamma) or \code{"B"} (beta / delta).}
#'   \item{region}{Character. Gene segment: \code{"V"}, \code{"D"}, or \code{"J"}.}
#'   \item{nucseq}{Character. Nucleotide sequence.}
#'   \item{alseq}{Character. Aligned protein sequence (gaps represented as \code{"."}).}
#'   \item{cdrs}{Character vector. CDR subsequences extracted from \code{alseq}.}
#'   \item{cdr_columns}{List of 2-element integer vectors. Start and end positions
#'     (1-indexed, inclusive) of each CDR in \code{alseq}.}
#'   \item{nucseq_offset}{Integer. 0-based reading frame offset (frame - 1).}
#'   \item{protseq}{Character. Protein sequence without gap characters.}
#'   \item{rep}{Character. Representative gene ID from exact loopseq neighbours (min by ID).}
#'   \item{mm1_rep}{Character. Representative gene ID from transitive mm1 loopseq neighbours.}
#'   \item{count_rep}{Character. Gene-level name (allele stripped); used for clone counting.}
#' }
#' @examples
#' \donttest{
#' # Load all organisms
#' all_g <- load_gene_database()
#' names(all_g)  # "human", "mouse", ...
#'
#' # Load one organism
#' human_genes <- load_gene_database("human")
#' human_genes[["TRAV1-1*01"]]$protseq
#' }
#' @seealso \code{\link{tcrdist_matrix}}, \code{\link{trim_allele_to_gene}}
#' @export
load_gene_database <- function(organism = NULL) {
    # Lazy-load and cache
    if (length(.tcrdistR_env$all_genes) == 0L) {
        .tcrdistR_env$all_genes <- .build_gene_database()
    }

    if (is.null(organism)) {
        return(.tcrdistR_env$all_genes)
    }

    if (!organism %in% names(.tcrdistR_env$all_genes)) {
        stop(sprintf(
            "Organism '%s' not found in gene database. Available: %s",
            organism,
            paste(sort(names(.tcrdistR_env$all_genes)), collapse = ", ")
        ))
    }

    .tcrdistR_env$all_genes[[organism]]
}
