# R wrapper for sparse TCRdist matrix computation.

# ---------------------------------------------------------------------------
# tcrdist_sparse
# ---------------------------------------------------------------------------

#' Compute sparse TCRdist matrix as a \code{dgCMatrix}
#'
#' Computes pairwise TCRdist distances for all pairs \code{(i, j)} where
#' \code{i < j} and the distance is at most \code{threshold}, then returns a
#' symmetric sparse matrix in compressed-column format (\code{dgCMatrix}).
#'
#' Three-stage early termination is used internally:
#' \enumerate{
#'   \item If the V-region distance sum alone exceeds \code{threshold}, skip.
#'   \item If V-region + CDR3-alpha distance exceeds \code{threshold}, skip.
#'   \item If the full distance exceeds \code{threshold}, skip.
#' }
#' For \code{threshold = Inf} all pairs are evaluated (equivalent to a full
#' dense matrix stored as sparse), which serves as a correctness check against
#' \code{\link{tcrdist_matrix}}.
#'
#' @param tcrs A \code{data.frame} with at least the following columns:
#'   \describe{
#'     \item{\code{va}}{Character. Alpha-chain V-gene allele.}
#'     \item{\code{cdr3a}}{Character. Alpha-chain CDR3 amino acid sequence.}
#'     \item{\code{vb}}{Character. Beta-chain V-gene allele.}
#'     \item{\code{cdr3b}}{Character. Beta-chain CDR3 amino acid sequence.}
#'   }
#' @param organism Character string. Organism key, e.g. \code{"human"}.
#' @param threshold Numeric. Maximum distance to include. Pairs with distance
#'   strictly greater than \code{threshold} are omitted. Must be >= 0.
#' @param weight_cdr3 Integer. CDR3 distance weight. Defaults to
#'   \code{WEIGHT_CDR3_REGION} (3L).
#' @param gap_penalty_cdr3 Integer. CDR3 gap penalty. Defaults to
#'   \code{GAP_PENALTY_CDR3_REGION} (12L).
#' @return A symmetric sparse matrix of class \code{dgCMatrix} with dimensions
#'   N x N.  Off-diagonal entries (i, j) and (j, i) are present for all pairs
#'   within the threshold.  The diagonal is structural zero (not stored).
#'   Returns an all-zero N x N \code{dgCMatrix} if no pairs satisfy the
#'   threshold.
#' @importFrom Matrix sparseMatrix
#' @examples
#' \donttest{
#' tcrs <- data.frame(
#'   va    = c("TRAV1-1*01", "TRAV1-1*01"),
#'   cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
#'   vb    = c("TRBV19*01", "TRBV19*01"),
#'   cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
#'   stringsAsFactors = FALSE
#' )
#' sp <- tcrdist_sparse(tcrs, "human", threshold = 50)
#' }
#' @export
tcrdist_sparse <- function(tcrs, organism, threshold,
                           weight_cdr3      = WEIGHT_CDR3_REGION,
                           gap_penalty_cdr3 = GAP_PENALTY_CDR3_REGION) {

    # ---- Input validation ---------------------------------------------------
    if (!is.data.frame(tcrs)) {
        stop("tcrdist_sparse: 'tcrs' must be a data.frame")
    }

    required_cols <- c("va", "cdr3a", "vb", "cdr3b")
    missing_cols  <- setdiff(required_cols, colnames(tcrs))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "tcrdist_sparse: missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
        threshold < 0) {
        stop("tcrdist_sparse: 'threshold' must be a non-negative numeric scalar")
    }

    n <- nrow(tcrs)
    if (n == 0L) {
        return(Matrix::sparseMatrix(i = integer(0L), j = integer(0L),
                                    x = numeric(0L),
                                    dims = c(0L, 0L)))
    }

    # ---- Coerce factors and check types ------------------------------------
    for (col in required_cols) {
        if (is.factor(tcrs[[col]])) tcrs[[col]] <- as.character(tcrs[[col]])
        if (!is.character(tcrs[[col]])) {
            stop(sprintf(
                "tcrdist_sparse: column '%s' must be character (or factor coercible to character)",
                col
            ))
        }
    }

    # ---- Check for NAs ------------------------------------------------------
    for (col in required_cols) {
        if (anyNA(tcrs[[col]])) {
            stop(sprintf(
                "tcrdist_sparse: column '%s' contains NA values", col
            ))
        }
    }

    # ---- Organism validation ------------------------------------------------
    if (!is.character(organism) || length(organism) != 1L || !nzchar(organism)) {
        stop("tcrdist_sparse: 'organism' must be a non-empty character string of length 1")
    }

    # ---- Build V-region distance matrices ----------------------------------
    v_dist_a <- .compute_v_region_distance_matrix(organism, "A")
    v_dist_b <- .compute_v_region_distance_matrix(organism, "B")

    # ---- Validate V genes --------------------------------------------------
    va_rownames <- rownames(v_dist_a)
    vb_rownames <- rownames(v_dist_b)

    unknown_va <- setdiff(unique(tcrs$va), va_rownames)
    if (length(unknown_va) > 0L) {
        stop(sprintf(
            "tcrdist_sparse: the following alpha V genes were not found for organism '%s': %s",
            organism, paste(unknown_va, collapse = ", ")
        ))
    }

    unknown_vb <- setdiff(unique(tcrs$vb), vb_rownames)
    if (length(unknown_vb) > 0L) {
        stop(sprintf(
            "tcrdist_sparse: the following beta V genes were not found for organism '%s': %s",
            organism, paste(unknown_vb, collapse = ", ")
        ))
    }

    # ---- Dispatch to C++ ---------------------------------------------------
    triplets <- rcpp_tcrdist_sparse(
        tcrs$va,
        tcrs$cdr3a,
        tcrs$vb,
        tcrs$cdr3b,
        v_dist_a,
        v_dist_b,
        threshold,
        as.integer(weight_cdr3),
        as.integer(gap_penalty_cdr3)
    )

    idx_names <- as.character(seq_len(n))

    # ---- Symmetrize COO triplets -------------------------------------------
    # C++ returns upper triangle only (i < j); mirror to get symmetric matrix.
    if (length(triplets$i) == 0L) {
        # No pairs within threshold: return zero sparse matrix
        return(Matrix::sparseMatrix(i = integer(0L), j = integer(0L),
                                    x = numeric(0L),
                                    dims    = c(n, n),
                                    dimnames = list(idx_names, idx_names)))
    }

    all_i <- c(triplets$i, triplets$j)
    all_j <- c(triplets$j, triplets$i)
    all_x <- c(triplets$x, triplets$x)

    Matrix::sparseMatrix(
        i        = all_i,
        j        = all_j,
        x        = all_x,
        dims     = c(n, n),
        dimnames = list(idx_names, idx_names)
    )
}
