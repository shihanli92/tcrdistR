# R wrapper for rectangular (query vs reference) TCRdist computation.

# ---------------------------------------------------------------------------
# tcrdist_rect
# ---------------------------------------------------------------------------

#' Compute rectangular TCRdist matrix between a query set and a reference set
#'
#' Computes an nq x nr matrix of paired-chain TCRdist distances between a
#' query set of \code{nq} TCRs and a reference set of \code{nr} TCRs.
#' Unlike \code{\link{tcrdist_matrix}}, the result is not symmetric because
#' query and reference need not be the same set.
#'
#' For each query-reference pair \code{(i, j)} the distance is:
#' \deqn{
#'   d(i,j) = v\_dist\_a[\mathrm{va\_query}_i, \mathrm{va\_ref}_j]
#'           + \mathrm{cdr3\_dist}(\mathrm{cdr3a\_query}_i, \mathrm{cdr3a\_ref}_j)
#'           + v\_dist\_b[\mathrm{vb\_query}_i, \mathrm{vb\_ref}_j]
#'           + \mathrm{cdr3\_dist}(\mathrm{cdr3b\_query}_i, \mathrm{cdr3b\_ref}_j)
#' }
#'
#' When \code{query} and \code{ref} are the same data.frame the result equals
#' \code{tcrdist_matrix(query, organism)} (no symmetry short-cut is taken so
#' the diagonal may differ by floating-point rounding from zero, but is
#' numerically identical when both data.frames are identical objects).
#'
#' @param query A \code{data.frame} with at least the following columns:
#'   \describe{
#'     \item{\code{va}}{Character. Alpha-chain V-gene allele, e.g. \code{"TRAV1-1*01"}.}
#'     \item{\code{cdr3a}}{Character. Alpha-chain CDR3 amino acid sequence.}
#'     \item{\code{vb}}{Character. Beta-chain V-gene allele, e.g. \code{"TRBV19*01"}.}
#'     \item{\code{cdr3b}}{Character. Beta-chain CDR3 amino acid sequence.}
#'   }
#' @param ref A \code{data.frame} with the same four columns as \code{query}.
#' @param organism Character string. Organism key understood by
#'   \code{load_gene_database}, e.g. \code{"human"} or \code{"mouse"}.
#' @param weight_cdr3 Integer. Weight applied to CDR3 distances. Defaults to
#'   \code{WEIGHT_CDR3_REGION} (3L).
#' @param gap_penalty_cdr3 Integer. Gap penalty for CDR3 alignments. Defaults
#'   to \code{GAP_PENALTY_CDR3_REGION} (12L).
#' @return A numeric matrix of dimensions nq x nr.  Row names are row indices
#'   of \code{query} as character strings; column names are row indices of
#'   \code{ref}.
#' @examples
#' \donttest{
#' tcrs <- data.frame(
#'   va    = c("TRAV1-1*01", "TRAV1-1*01"),
#'   cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
#'   vb    = c("TRBV19*01", "TRBV19*01"),
#'   cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
#'   stringsAsFactors = FALSE
#' )
#' # Self-comparison == tcrdist_matrix
#' rect <- tcrdist_rect(tcrs, tcrs, "human")
#' mat  <- tcrdist_matrix(tcrs, "human")
#' }
#' @seealso \code{\link{tcrdist_matrix}}, \code{\link{tcrdist_sparse}}, \code{\link{tcrdist_join}}
#' @export
tcrdist_rect <- function(query, ref, organism,
                         weight_cdr3      = WEIGHT_CDR3_REGION,
                         gap_penalty_cdr3 = GAP_PENALTY_CDR3_REGION) {

    # ---- Input validation: query --------------------------------------------
    if (!is.data.frame(query)) {
        stop("tcrdist_rect: 'query' must be a data.frame")
    }

    required_cols <- c("va", "cdr3a", "vb", "cdr3b")
    missing_query <- setdiff(required_cols, colnames(query))
    if (length(missing_query) > 0L) {
        stop(sprintf(
            "tcrdist_rect: 'query' is missing required columns: %s",
            paste(missing_query, collapse = ", ")
        ))
    }

    # ---- Input validation: ref ----------------------------------------------
    if (!is.data.frame(ref)) {
        stop("tcrdist_rect: 'ref' must be a data.frame")
    }

    missing_ref <- setdiff(required_cols, colnames(ref))
    if (length(missing_ref) > 0L) {
        stop(sprintf(
            "tcrdist_rect: 'ref' is missing required columns: %s",
            paste(missing_ref, collapse = ", ")
        ))
    }

    nq <- nrow(query)
    nr <- nrow(ref)

    if (nq == 0L || nr == 0L) {
        return(matrix(numeric(0L), nrow = nq, ncol = nr))
    }

    # ---- Coerce factors and check types ------------------------------------
    for (col in required_cols) {
        if (is.factor(query[[col]])) query[[col]] <- as.character(query[[col]])
        if (!is.character(query[[col]])) {
            stop(sprintf(
                "tcrdist_rect: query column '%s' must be character (or factor coercible to character)",
                col
            ))
        }
        if (is.factor(ref[[col]])) ref[[col]] <- as.character(ref[[col]])
        if (!is.character(ref[[col]])) {
            stop(sprintf(
                "tcrdist_rect: ref column '%s' must be character (or factor coercible to character)",
                col
            ))
        }
    }

    # ---- Check for NAs ------------------------------------------------------
    for (col in required_cols) {
        if (anyNA(query[[col]])) {
            stop(sprintf(
                "tcrdist_rect: query column '%s' contains NA values", col
            ))
        }
        if (anyNA(ref[[col]])) {
            stop(sprintf(
                "tcrdist_rect: ref column '%s' contains NA values", col
            ))
        }
    }

    # ---- Organism validation ------------------------------------------------
    if (!is.character(organism) || length(organism) != 1L || !nzchar(organism)) {
        stop("tcrdist_rect: 'organism' must be a non-empty character string of length 1")
    }

    # ---- Build V-region distance matrices ----------------------------------
    v_dist_a <- .compute_v_region_distance_matrix(organism, "A")
    v_dist_b <- .compute_v_region_distance_matrix(organism, "B")

    # ---- Validate V genes exist in both query and ref ----------------------
    va_rownames <- rownames(v_dist_a)
    vb_rownames <- rownames(v_dist_b)

    unknown_q_va <- setdiff(unique(query$va), va_rownames)
    if (length(unknown_q_va) > 0L) {
        stop(sprintf(
            "tcrdist_rect: the following query alpha V genes were not found for organism '%s': %s",
            organism, paste(unknown_q_va, collapse = ", ")
        ))
    }

    unknown_r_va <- setdiff(unique(ref$va), va_rownames)
    if (length(unknown_r_va) > 0L) {
        stop(sprintf(
            "tcrdist_rect: the following ref alpha V genes were not found for organism '%s': %s",
            organism, paste(unknown_r_va, collapse = ", ")
        ))
    }

    unknown_q_vb <- setdiff(unique(query$vb), vb_rownames)
    if (length(unknown_q_vb) > 0L) {
        stop(sprintf(
            "tcrdist_rect: the following query beta V genes were not found for organism '%s': %s",
            organism, paste(unknown_q_vb, collapse = ", ")
        ))
    }

    unknown_r_vb <- setdiff(unique(ref$vb), vb_rownames)
    if (length(unknown_r_vb) > 0L) {
        stop(sprintf(
            "tcrdist_rect: the following ref beta V genes were not found for organism '%s': %s",
            organism, paste(unknown_r_vb, collapse = ", ")
        ))
    }

    # ---- Dispatch to C++ ---------------------------------------------------
    dist_mat <- rcpp_tcrdist_rect(
        query$va,
        query$cdr3a,
        query$vb,
        query$cdr3b,
        ref$va,
        ref$cdr3a,
        ref$vb,
        ref$cdr3b,
        v_dist_a,
        v_dist_b,
        as.integer(weight_cdr3),
        as.integer(gap_penalty_cdr3)
    )

    # ---- Set row/column names ----------------------------------------------
    rownames(dist_mat) <- as.character(seq_len(nq))
    colnames(dist_mat) <- as.character(seq_len(nr))

    dist_mat
}
