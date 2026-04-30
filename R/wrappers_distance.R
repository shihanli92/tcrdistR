# Thin R wrappers that validate inputs and dispatch to C++ rcpp_* functions.


# ---------------------------------------------------------------------------
# .resolve_components  (internal)
# ---------------------------------------------------------------------------

#' Resolve a components specification into per-component flags
#'
#' @param components Character.  A single preset string or a vector of
#'   individual component names.
#' @return A named logical list with elements \code{va}, \code{cdr3a},
#'   \code{vb}, \code{cdr3b}.
#' @keywords internal
.resolve_components <- function(components) {
    valid_terms <- c("va", "cdr3a", "vb", "cdr3b")
    presets <- list(
        all      = valid_terms,
        cdr3     = c("cdr3a", "cdr3b"),
        v_region = c("va", "vb"),
        alpha    = c("va", "cdr3a"),
        beta     = c("vb", "cdr3b")
    )

    if (length(components) == 1L && components %in% names(presets)) {
        active <- presets[[components]]
    } else {
        unknown <- setdiff(components, valid_terms)
        if (length(unknown) > 0L) {
            stop(sprintf(
                "Unknown components: %s. Valid terms: %s. Presets: %s",
                paste(unknown, collapse = ", "),
                paste(valid_terms, collapse = ", "),
                paste(names(presets), collapse = ", ")
            ), call. = FALSE)
        }
        active <- components
    }

    list(
        va    = "va"    %in% active,
        cdr3a = "cdr3a" %in% active,
        vb    = "vb"    %in% active,
        cdr3b = "cdr3b" %in% active
    )
}


#' Fill missing chain columns for single-chain input
#'
#' Detects whether both chains are present and, if not, fills the missing
#' chain with dummy values so that downstream code can proceed unchanged.
#' The dummy V-gene is the first V-gene for that chain in the gene database;
#' the dummy CDR3 is a minimal valid sequence.  When combined with the
#' \code{components} parameter (auto-set to \code{"alpha"} or \code{"beta"}),
#' the dummy chain contributes exactly zero to the distance.
#'
#' @param tcrs A \code{data.frame}.
#' @param organism Character string passed to \code{load_gene_database}.
#' @return A list with elements \code{tcrs} (possibly augmented) and
#'   \code{chain} (\code{"AB"}, \code{"A"}, or \code{"B"}).
#' @keywords internal
.fill_missing_chain <- function(tcrs, organism) {
    has_alpha <- all(c("va", "cdr3a") %in% colnames(tcrs))
    has_beta  <- all(c("vb", "cdr3b") %in% colnames(tcrs))

    if (has_alpha && has_beta) {
        return(list(tcrs = tcrs, chain = "AB"))
    }
    if (!has_alpha && !has_beta) {
        stop(
            "tcrs must contain at least alpha (va, cdr3a) or beta (vb, cdr3b) columns",
            call. = FALSE
        )
    }

    all_genes_org <- load_gene_database(organism)
    ids <- names(all_genes_org)

    if (!has_alpha) {
        # Beta-only: fill alpha with dummies
        dummy_va <- ids[vapply(ids, function(id) {
            g <- all_genes_org[[id]]
            g$chain == "A" && g$region == "V"
        }, logical(1L))][1L]
        tcrs$va    <- dummy_va
        tcrs$cdr3a <- "CAAAAF"
        return(list(tcrs = tcrs, chain = "B"))
    }

    # Alpha-only: fill beta with dummies
    dummy_vb <- ids[vapply(ids, function(id) {
        g <- all_genes_org[[id]]
        g$chain == "B" && g$region == "V"
    }, logical(1L))][1L]
    tcrs$vb    <- dummy_vb
    tcrs$cdr3b <- "CASSF"
    list(tcrs = tcrs, chain = "A")
}


#' Create a zero V-region distance matrix with matching gene names
#'
#' @param v_dist A named square \code{NumericMatrix} from
#'   \code{.compute_v_region_distance_matrix()}.
#' @return A zero matrix with the same dimensions and row/column names.
#' @keywords internal
.zero_v_dist_matrix <- function(v_dist) {
    n <- nrow(v_dist)
    m <- matrix(0, nrow = n, ncol = n)
    rownames(m) <- rownames(v_dist)
    colnames(m) <- colnames(v_dist)
    m
}


# ---------------------------------------------------------------------------
# weighted_cdr3_distance
# ---------------------------------------------------------------------------

#' Weighted CDR3 distance between two CDR3 sequences
#'
#' Computes the TCRdist CDR3 component distance between two amino acid CDR3
#' sequences. Accounts for length differences via gap penalties and uses BSD4
#' substitution scores on the aligned positions. Dispatches to a C++
#' implementation for performance.
#'
#' @param seq1 Character string of length 1. First CDR3 amino acid sequence.
#'   Must be non-empty.
#' @param seq2 Character string of length 1. Second CDR3 amino acid sequence.
#'   Must be non-empty.
#' @param weight Integer. Weight applied to the alignment score component.
#'   Defaults to \code{WEIGHT_CDR3_REGION} (3L).
#' @param gap_penalty Integer. Penalty per gap character in length differences.
#'   Defaults to \code{GAP_PENALTY_CDR3_REGION} (12L).
#' @return A numeric scalar: the weighted CDR3 distance.
#' @examples
#' \donttest{
#' weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSYEQYF")
#' weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSSYEQYF")  # 0
#' }
#' @seealso \code{\link{tcrdist_matrix}}, \code{\link{bsd4_matrix}}
#' @export
weighted_cdr3_distance <- function(seq1, seq2,
                                   weight      = WEIGHT_CDR3_REGION,
                                   gap_penalty = GAP_PENALTY_CDR3_REGION) {
    if (!is.character(seq1) || length(seq1) != 1L || is.na(seq1) || !nzchar(seq1)) {
        stop("weighted_cdr3_distance: 'seq1' must be a non-empty, non-NA character string")
    }
    if (!is.character(seq2) || length(seq2) != 1L || is.na(seq2) || !nzchar(seq2)) {
        stop("weighted_cdr3_distance: 'seq2' must be a non-empty, non-NA character string")
    }

    rcpp_weighted_cdr3_distance(seq1, seq2,
                                as.integer(weight),
                                as.integer(gap_penalty))
}


# ---------------------------------------------------------------------------
# tcrdist_matrix
# ---------------------------------------------------------------------------

#' Compute pairwise TCRdist distance matrix
#'
#' Computes the full N x N symmetric matrix of paired-chain TCRdist distances
#' for a collection of TCRs. Each off-diagonal entry \code{[i, j]} equals the
#' sum of V-alpha, CDR3-alpha, V-beta, and CDR3-beta component distances
#' between TCR \code{i} and TCR \code{j}. Computation is dispatched to a C++
#' implementation for performance.
#'
#' The diagonal is zero (distance of a TCR to itself). The matrix is symmetric
#' by construction.
#'
#' @param tcrs A \code{data.frame} with columns for one or both TCR chains.
#'   For paired alpha-beta input, requires \code{va}, \code{cdr3a}, \code{vb},
#'   \code{cdr3b}.  For single-chain input, only the columns for one chain are
#'   needed (e.g., \code{vb} and \code{cdr3b} for beta-only).  The missing
#'   chain is filled with dummy values internally, and \code{components} is
#'   automatically set to \code{"alpha"} or \code{"beta"}.
#' @param organism Character string. Organism key understood by
#'   \code{load_gene_database}, e.g. \code{"human"} or \code{"mouse"}.
#' @param components Character.  Which distance components to include.
#'   Presets: \code{"all"} (default), \code{"cdr3"} (CDR3 only),
#'   \code{"v_region"} (CDR1+CDR2+CDR2.5 only), \code{"alpha"} (alpha chain),
#'   \code{"beta"} (beta chain).  Or a character vector of individual terms:
#'   \code{"va"}, \code{"cdr3a"}, \code{"vb"}, \code{"cdr3b"}.
#'   For single-chain input, defaults to the present chain.
#' @param weight_cdr3 Integer. Weight applied to CDR3 distances. Defaults to
#'   \code{WEIGHT_CDR3_REGION} (3L).
#' @param gap_penalty_cdr3 Integer. Gap penalty for CDR3 alignments. Defaults
#'   to \code{GAP_PENALTY_CDR3_REGION} (12L).
#' @return A symmetric numeric matrix of dimensions N x N where N is
#'   \code{nrow(tcrs)}. Row and column names are the row indices of \code{tcrs}
#'   as character strings.
#' @examples
#' \donttest{
#' tcrs <- data.frame(
#'   va    = c("TRAV1-1*01", "TRAV1-1*01"),
#'   cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
#'   vb    = c("TRBV19*01", "TRBV19*01"),
#'   cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
#'   stringsAsFactors = FALSE
#' )
#' mat <- tcrdist_matrix(tcrs, "human")
#'
#' # CDR3-only distance
#' mat_cdr3 <- tcrdist_matrix(tcrs, "human", components = "cdr3")
#'
#' # Beta-only (single-chain) input
#' beta_only <- tcrs[, c("vb", "cdr3b")]
#' mat_beta <- tcrdist_matrix(beta_only, "human")
#' }
#' @seealso \code{\link{tcrdist_sparse}}, \code{\link{tcrdist_rect}}, \code{\link{tcrdist_knn}}, \code{\link{TCRrep}}
#' @export
tcrdist_matrix <- function(tcrs, organism,
                           components       = "all",
                           weight_cdr3      = WEIGHT_CDR3_REGION,
                           gap_penalty_cdr3 = GAP_PENALTY_CDR3_REGION) {
    # ---- Input validation ---------------------------------------------------
    if (!is.data.frame(tcrs)) {
        stop("tcrdist_matrix: 'tcrs' must be a data.frame")
    }

    # ---- Single-chain detection -----------------------------------------------
    filled <- .fill_missing_chain(tcrs, organism)
    tcrs <- filled$tcrs
    if (filled$chain != "AB" && identical(components, "all")) {
        components <- if (filled$chain == "A") "alpha" else "beta"
    }

    required_cols <- c("va", "cdr3a", "vb", "cdr3b")

    n <- nrow(tcrs)
    if (n == 0L) {
        return(matrix(numeric(0L), nrow = 0L, ncol = 0L))
    }

    if (n > 20000L) {
        warning(sprintf(
            "tcrdist_matrix: N=%d requires a full %d x %d distance matrix (%.1f GB). Consider tcrdist_sparse() or tcrdist_radius_neighbors() for large datasets.",
            n, n, n, as.double(n) * n * 8 / 1e9
        ))
    }

    # Ensure character columns (not factors)
    for (col in required_cols) {
        if (is.factor(tcrs[[col]])) {
            tcrs[[col]] <- as.character(tcrs[[col]])
        }
        if (!is.character(tcrs[[col]])) {
            stop(sprintf(
                "tcrdist_matrix: column '%s' must be character (or factor coercible to character)",
                col
            ))
        }
    }

    # Check for NAs
    for (col in required_cols) {
        if (anyNA(tcrs[[col]])) {
            stop(sprintf(
                "tcrdist_matrix: column '%s' contains NA values",
                col
            ))
        }
    }

    if (!is.character(organism) || length(organism) != 1L || !nzchar(organism)) {
        stop("tcrdist_matrix: 'organism' must be a non-empty character string of length 1")
    }

    # ---- Resolve components -------------------------------------------------
    comp <- .resolve_components(components)

    # ---- Build V-region distance matrices via C++ --------------------------
    v_dist_a <- .compute_v_region_distance_matrix(organism, "A")
    v_dist_b <- .compute_v_region_distance_matrix(organism, "B")

    # ---- Validate that all V genes exist in the distance matrices ----------
    va_rownames <- rownames(v_dist_a)
    vb_rownames <- rownames(v_dist_b)

    unknown_va <- setdiff(unique(tcrs$va), va_rownames)
    if (length(unknown_va) > 0L) {
        stop(sprintf(
            "tcrdist_matrix: the following alpha V genes were not found in the distance matrix for organism '%s': %s",
            organism,
            paste(unknown_va, collapse = ", ")
        ))
    }

    unknown_vb <- setdiff(unique(tcrs$vb), vb_rownames)
    if (length(unknown_vb) > 0L) {
        stop(sprintf(
            "tcrdist_matrix: the following beta V genes were not found in the distance matrix for organism '%s': %s",
            organism,
            paste(unknown_vb, collapse = ", ")
        ))
    }

    # ---- Apply component selection -----------------------------------------
    if (!comp$va)    v_dist_a <- .zero_v_dist_matrix(v_dist_a)
    if (!comp$vb)    v_dist_b <- .zero_v_dist_matrix(v_dist_b)

    w_a  <- if (comp$cdr3a) as.integer(weight_cdr3)      else 0L
    gp_a <- if (comp$cdr3a) as.integer(gap_penalty_cdr3) else 0L
    w_b  <- if (comp$cdr3b) as.integer(weight_cdr3)      else 0L
    gp_b <- if (comp$cdr3b) as.integer(gap_penalty_cdr3) else 0L

    # ---- Dispatch to C++ ---------------------------------------------------
    dist_mat <- rcpp_tcrdist_matrix(
        tcrs$va,
        tcrs$cdr3a,
        tcrs$vb,
        tcrs$cdr3b,
        v_dist_a,
        v_dist_b,
        w_a, gp_a,
        w_b, gp_b
    )

    # ---- Set row/column names ----------------------------------------------
    idx_names           <- as.character(seq_len(n))
    rownames(dist_mat)  <- idx_names
    colnames(dist_mat)  <- idx_names

    dist_mat
}


# ---------------------------------------------------------------------------
# bsd4_matrix
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# hamming_distance
# ---------------------------------------------------------------------------

#' Hamming distance between two CDR3 sequences
#'
#' Counts the number of positional mismatches between two equal-length amino
#' acid sequences. Returns \code{-1L} if the sequences have different lengths.
#'
#' @param a Character string. First CDR3 amino acid sequence.
#' @param b Character string. Second CDR3 amino acid sequence.
#' @return An integer: the number of mismatches, or \code{-1L} if lengths
#'   differ.
#' @examples
#' hamming_distance("CASSI", "CASSK")  # 1
#' hamming_distance("CASSI", "CASSI")  # 0
#' hamming_distance("CASSI", "CASSILY")  # -1 (different lengths)
#' @seealso \code{\link{hamming_matrix}}, \code{\link{weighted_cdr3_distance}}
#' @export
hamming_distance <- function(a, b) {
    stopifnot(
        is.character(a), length(a) == 1L, !is.na(a),
        is.character(b), length(b) == 1L, !is.na(b)
    )
    rcpp_hamming_distance(a, b)
}


# ---------------------------------------------------------------------------
# hamming_matrix
# ---------------------------------------------------------------------------

#' Pairwise Hamming distance matrix for CDR3 sequences
#'
#' Computes an N x N integer matrix of Hamming distances between CDR3
#' sequences. For pairs of equal length, the distance is the number of
#' mismatches. For pairs of unequal length, the distance is set to the
#' length of the longer sequence (maximum penalty).
#'
#' @param cdr3_seqs Character vector. CDR3 amino acid sequences.
#' @return An integer matrix of dimensions N x N.
#' @examples
#' seqs <- c("CASSI", "CASSK", "CASRL")
#' hamming_matrix(seqs)
#' @seealso \code{\link{hamming_distance}}, \code{\link{tcrdist_matrix}}
#' @export
hamming_matrix <- function(cdr3_seqs) {
    stopifnot(is.character(cdr3_seqs), length(cdr3_seqs) >= 1L)
    if (anyNA(cdr3_seqs)) {
        stop("hamming_matrix: cdr3_seqs must not contain NA values")
    }
    rcpp_hamming_matrix(cdr3_seqs)
}


# ---------------------------------------------------------------------------
# bsd4_matrix
# ---------------------------------------------------------------------------

#' Retrieve the BSD4 substitution matrix
#'
#' Returns the BSD4 (BLOSUM-derived substitution distance 4) matrix used
#' internally for amino acid distance scoring in TCRdist calculations. The
#' matrix is built in C++ and returned as a named numeric matrix in R.
#'
#' @return A named numeric matrix of amino acid substitution distances. Row
#'   and column names are the 20 standard single-letter amino acid codes in
#'   \code{AMINO_ACIDS} order.
#' @examples
#' \donttest{
#' bsd4 <- bsd4_matrix()
#' bsd4["A", "A"]  # 0
#' bsd4["A", "G"]  # small positive value
#' }
#' @seealso \code{\link{weighted_cdr3_distance}}, \code{\link{AMINO_ACIDS}}
#' @export
bsd4_matrix <- function() {
    rcpp_build_bsd4()
}
