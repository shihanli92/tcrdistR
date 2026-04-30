# R wrappers for KNN and radius neighbor functions.

# ---------------------------------------------------------------------------
# Internal helper: default group vector (no masking)
# ---------------------------------------------------------------------------
# When agroups/bgroups is NULL, assign each TCR its own unique group so that
# same-group masking has no effect.

.default_groups <- function(n) seq_len(n)

# ---------------------------------------------------------------------------
# tcrdist_knn
# ---------------------------------------------------------------------------

#' K-nearest-neighbors by TCRdist with optional group masking
#'
#' For each of the N input TCRs, finds the K nearest neighbors by TCRdist
#' distance.  TCRs sharing the same \code{agroups} or \code{bgroups} value are
#' excluded from each other's neighborhood (same-group masking).  When
#' \code{agroups} and \code{bgroups} are \code{NULL} (default), every TCR has
#' its own unique group so no masking occurs.
#'
#' @param tcrs A \code{data.frame} with at least the following columns:
#'   \describe{
#'     \item{\code{va}}{Character. Alpha-chain V-gene allele.}
#'     \item{\code{cdr3a}}{Character. Alpha-chain CDR3 amino acid sequence.}
#'     \item{\code{vb}}{Character. Beta-chain V-gene allele.}
#'     \item{\code{cdr3b}}{Character. Beta-chain CDR3 amino acid sequence.}
#'   }
#' @param organism Character string. Organism key, e.g. \code{"human"}.
#' @param K Integer. Number of nearest neighbors to return per TCR.
#'   Must satisfy \code{1 <= K <= N - 1}.
#' @param agroups Integer vector of length N, or \code{NULL}. Alpha-chain
#'   group assignments.  TCRs with the same value are masked from each other.
#'   \code{NULL} assigns each TCR a unique group (no masking).
#' @param bgroups Integer vector of length N, or \code{NULL}. Beta-chain
#'   group assignments. Same semantics as \code{agroups}.
#' @param sort_nbrs Logical. If \code{TRUE} (default), sort each row's K
#'   neighbors by ascending distance.
#' @param components Character.  Which distance components to include.
#'   See \code{\link{tcrdist_matrix}} for details.  Default \code{"all"}.
#' @param weight_cdr3 Integer. CDR3 distance weight. Defaults to
#'   \code{WEIGHT_CDR3_REGION} (3L).
#' @param gap_penalty_cdr3 Integer. CDR3 gap penalty. Defaults to
#'   \code{GAP_PENALTY_CDR3_REGION} (12L).
#' @return A \code{list} with two elements:
#'   \describe{
#'     \item{\code{knn_indices}}{Integer matrix (N x K). 1-based row indices
#'       of the K nearest neighbors for each TCR.}
#'     \item{\code{knn_distances}}{Numeric matrix (N x K). Corresponding
#'       TCRdist distances.}
#'   }
#' @examples
#' \donttest{
#' tcrs <- data.frame(
#'   va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
#'   cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
#'   vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
#'   cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
#'   stringsAsFactors = FALSE
#' )
#' knn <- tcrdist_knn(tcrs, "human", K = 1L)
#' }
#' @seealso \code{\link{tcrdist_radius_neighbors}}, \code{\link{knn_from_matrix}}, \code{\link{knn_from_pca}}
#' @export
tcrdist_knn <- function(tcrs, organism, K,
                        agroups          = NULL,
                        bgroups          = NULL,
                        sort_nbrs        = TRUE,
                        components       = "all",
                        weight_cdr3      = WEIGHT_CDR3_REGION,
                        gap_penalty_cdr3 = GAP_PENALTY_CDR3_REGION) {

    # ---- Input validation ---------------------------------------------------
    if (!is.data.frame(tcrs)) {
        stop("tcrdist_knn: 'tcrs' must be a data.frame")
    }

    required_cols <- c("va", "cdr3a", "vb", "cdr3b")
    missing_cols  <- setdiff(required_cols, colnames(tcrs))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "tcrdist_knn: missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    n <- nrow(tcrs)
    if (n == 0L) {
        stop("tcrdist_knn: 'tcrs' is empty")
    }

    K <- as.integer(K)
    if (length(K) != 1L || is.na(K) || K < 1L || K >= n) {
        stop(sprintf(
            "tcrdist_knn: K must be an integer in [1, N-1]; got K=%d, N=%d", K, n
        ))
    }

    # ---- Default groups (unique per TCR = no masking) ----------------------
    if (is.null(agroups)) agroups <- .default_groups(n)
    if (is.null(bgroups)) bgroups <- .default_groups(n)

    agroups <- as.integer(agroups)
    bgroups <- as.integer(bgroups)

    if (length(agroups) != n) {
        stop(sprintf(
            "tcrdist_knn: agroups length (%d) != nrow(tcrs) (%d)",
            length(agroups), n
        ))
    }
    if (length(bgroups) != n) {
        stop(sprintf(
            "tcrdist_knn: bgroups length (%d) != nrow(tcrs) (%d)",
            length(bgroups), n
        ))
    }

    # ---- Coerce factors and check types ------------------------------------
    for (col in required_cols) {
        if (is.factor(tcrs[[col]])) tcrs[[col]] <- as.character(tcrs[[col]])
        if (!is.character(tcrs[[col]])) {
            stop(sprintf(
                "tcrdist_knn: column '%s' must be character (or factor coercible to character)",
                col
            ))
        }
    }

    # ---- Check for NAs ------------------------------------------------------
    for (col in required_cols) {
        if (anyNA(tcrs[[col]])) {
            stop(sprintf("tcrdist_knn: column '%s' contains NA values", col))
        }
    }

    # ---- Organism validation ------------------------------------------------
    if (!is.character(organism) || length(organism) != 1L || !nzchar(organism)) {
        stop("tcrdist_knn: 'organism' must be a non-empty character string of length 1")
    }

    # ---- Resolve components -------------------------------------------------
    comp <- .resolve_components(components)

    # ---- Build V-region distance matrices ----------------------------------
    v_dist_a <- .compute_v_region_distance_matrix(organism, "A")
    v_dist_b <- .compute_v_region_distance_matrix(organism, "B")

    # ---- Validate V genes --------------------------------------------------
    unknown_va <- setdiff(unique(tcrs$va), rownames(v_dist_a))
    if (length(unknown_va) > 0L) {
        stop(sprintf(
            "tcrdist_knn: the following alpha V genes were not found for organism '%s': %s",
            organism, paste(unknown_va, collapse = ", ")
        ))
    }
    unknown_vb <- setdiff(unique(tcrs$vb), rownames(v_dist_b))
    if (length(unknown_vb) > 0L) {
        stop(sprintf(
            "tcrdist_knn: the following beta V genes were not found for organism '%s': %s",
            organism, paste(unknown_vb, collapse = ", ")
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
    rcpp_tcrdist_knn(
        tcrs$va,
        tcrs$cdr3a,
        tcrs$vb,
        tcrs$cdr3b,
        v_dist_a,
        v_dist_b,
        K,
        agroups,
        bgroups,
        isTRUE(sort_nbrs),
        w_a, gp_a,
        w_b, gp_b
    )
}


# ---------------------------------------------------------------------------
# tcrdist_radius_neighbors
# ---------------------------------------------------------------------------

#' Radius-based TCRdist neighbor search with optional group masking
#'
#' For each of the N input TCRs, finds all other TCRs within \code{radius}
#' TCRdist distance.  TCRs sharing the same \code{agroups} or \code{bgroups}
#' value are excluded (same-group masking).
#'
#' @param tcrs A \code{data.frame} with at least the following columns:
#'   \describe{
#'     \item{\code{va}}{Character. Alpha-chain V-gene allele.}
#'     \item{\code{cdr3a}}{Character. Alpha-chain CDR3 amino acid sequence.}
#'     \item{\code{vb}}{Character. Beta-chain V-gene allele.}
#'     \item{\code{cdr3b}}{Character. Beta-chain CDR3 amino acid sequence.}
#'   }
#' @param organism Character string. Organism key, e.g. \code{"human"}.
#' @param radius Numeric. Search radius (inclusive). Pairs with distance <=
#'   \code{radius} are returned. Must be >= 0.
#' @param agroups Integer vector of length N, or \code{NULL}. Alpha-chain
#'   group assignments.  \code{NULL} assigns each TCR a unique group (no masking).
#' @param bgroups Integer vector of length N, or \code{NULL}. Beta-chain
#'   group assignments.
#' @param components Character.  Which distance components to include.
#'   See \code{\link{tcrdist_matrix}} for details.  Default \code{"all"}.
#' @param weight_cdr3 Integer. CDR3 distance weight. Defaults to
#'   \code{WEIGHT_CDR3_REGION} (3L).
#' @param gap_penalty_cdr3 Integer. CDR3 gap penalty. Defaults to
#'   \code{GAP_PENALTY_CDR3_REGION} (12L).
#' @return A \code{list} of length N.  Each element \code{[[i]]} is a
#'   \code{list} with:
#'   \describe{
#'     \item{\code{indices}}{Integer vector. 1-based indices of neighbors
#'       within \code{radius}.}
#'     \item{\code{distances}}{Numeric vector. Corresponding distances.}
#'   }
#' @examples
#' \donttest{
#' tcrs <- data.frame(
#'   va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
#'   cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
#'   vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
#'   cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
#'   stringsAsFactors = FALSE
#' )
#' nbrs <- tcrdist_radius_neighbors(tcrs, "human", radius = 50)
#' }
#' @seealso \code{\link{tcrdist_knn}}, \code{\link{find_clumping}}
#' @export
tcrdist_radius_neighbors <- function(tcrs, organism, radius,
                                     agroups          = NULL,
                                     bgroups          = NULL,
                                     components       = "all",
                                     weight_cdr3      = WEIGHT_CDR3_REGION,
                                     gap_penalty_cdr3 = GAP_PENALTY_CDR3_REGION) {

    # ---- Input validation ---------------------------------------------------
    if (!is.data.frame(tcrs)) {
        stop("tcrdist_radius_neighbors: 'tcrs' must be a data.frame")
    }

    required_cols <- c("va", "cdr3a", "vb", "cdr3b")
    missing_cols  <- setdiff(required_cols, colnames(tcrs))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "tcrdist_radius_neighbors: missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) || radius < 0) {
        stop("tcrdist_radius_neighbors: 'radius' must be a non-negative numeric scalar")
    }

    n <- nrow(tcrs)
    if (n == 0L) {
        return(list())
    }

    # ---- Default groups (unique per TCR = no masking) ----------------------
    if (is.null(agroups)) agroups <- .default_groups(n)
    if (is.null(bgroups)) bgroups <- .default_groups(n)

    agroups <- as.integer(agroups)
    bgroups <- as.integer(bgroups)

    if (length(agroups) != n) {
        stop(sprintf(
            "tcrdist_radius_neighbors: agroups length (%d) != nrow(tcrs) (%d)",
            length(agroups), n
        ))
    }
    if (length(bgroups) != n) {
        stop(sprintf(
            "tcrdist_radius_neighbors: bgroups length (%d) != nrow(tcrs) (%d)",
            length(bgroups), n
        ))
    }

    # ---- Coerce factors and check types ------------------------------------
    for (col in required_cols) {
        if (is.factor(tcrs[[col]])) tcrs[[col]] <- as.character(tcrs[[col]])
        if (!is.character(tcrs[[col]])) {
            stop(sprintf(
                "tcrdist_radius_neighbors: column '%s' must be character (or factor coercible to character)",
                col
            ))
        }
    }

    # ---- Check for NAs ------------------------------------------------------
    for (col in required_cols) {
        if (anyNA(tcrs[[col]])) {
            stop(sprintf("tcrdist_radius_neighbors: column '%s' contains NA values", col))
        }
    }

    # ---- Organism validation ------------------------------------------------
    if (!is.character(organism) || length(organism) != 1L || !nzchar(organism)) {
        stop("tcrdist_radius_neighbors: 'organism' must be a non-empty character string of length 1")
    }

    # ---- Resolve components -------------------------------------------------
    comp <- .resolve_components(components)

    # ---- Build V-region distance matrices ----------------------------------
    v_dist_a <- .compute_v_region_distance_matrix(organism, "A")
    v_dist_b <- .compute_v_region_distance_matrix(organism, "B")

    # ---- Validate V genes --------------------------------------------------
    unknown_va <- setdiff(unique(tcrs$va), rownames(v_dist_a))
    if (length(unknown_va) > 0L) {
        stop(sprintf(
            "tcrdist_radius_neighbors: the following alpha V genes were not found for organism '%s': %s",
            organism, paste(unknown_va, collapse = ", ")
        ))
    }
    unknown_vb <- setdiff(unique(tcrs$vb), rownames(v_dist_b))
    if (length(unknown_vb) > 0L) {
        stop(sprintf(
            "tcrdist_radius_neighbors: the following beta V genes were not found for organism '%s': %s",
            organism, paste(unknown_vb, collapse = ", ")
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
    rcpp_tcrdist_radius_neighbors(
        tcrs$va,
        tcrs$cdr3a,
        tcrs$vb,
        tcrs$cdr3b,
        v_dist_a,
        v_dist_b,
        radius,
        agroups,
        bgroups,
        w_a, gp_a,
        w_b, gp_b
    )
}


# ---------------------------------------------------------------------------
# knn_from_matrix
# ---------------------------------------------------------------------------

#' K-nearest-neighbors from a precomputed distance matrix
#'
#' Extracts K nearest neighbors for each row of a precomputed N x N distance
#' matrix \code{D}, with optional same-group masking.
#'
#' @param D Numeric matrix (N x N). Precomputed pairwise distance matrix.
#'   Must be square.
#' @param K Integer. Number of nearest neighbors to return.
#'   Must satisfy \code{1 <= K <= N - 1}.
#' @param agroups Integer vector of length N, or \code{NULL}. Alpha-chain
#'   group assignments. \code{NULL} assigns unique groups (no masking).
#' @param bgroups Integer vector of length N, or \code{NULL}. Beta-chain
#'   group assignments.
#' @param sort_nbrs Logical. If \code{TRUE} (default), sort K neighbors by
#'   ascending distance.
#' @return A \code{list} with two elements:
#'   \describe{
#'     \item{\code{knn_indices}}{Integer matrix (N x K). 1-based indices.}
#'     \item{\code{knn_distances}}{Numeric matrix (N x K). Distances.}
#'   }
#' @examples
#' \donttest{
#' tcrs <- data.frame(
#'   va    = c("TRAV1-1*01", "TRAV1-1*01"),
#'   cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
#'   vb    = c("TRBV19*01", "TRBV19*01"),
#'   cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
#'   stringsAsFactors = FALSE
#' )
#' D   <- tcrdist_matrix(tcrs, "human")
#' knn <- knn_from_matrix(D, K = 1L)
#' }
#' @seealso \code{\link{knn_from_pca}}, \code{\link{tcrdist_knn}}
#' @export
knn_from_matrix <- function(D, K,
                             agroups   = NULL,
                             bgroups   = NULL,
                             sort_nbrs = TRUE) {

    # ---- Input validation ---------------------------------------------------
    if (!is.matrix(D) || !is.numeric(D)) {
        stop("knn_from_matrix: 'D' must be a numeric matrix")
    }
    n <- nrow(D)
    if (ncol(D) != n) {
        stop(sprintf(
            "knn_from_matrix: 'D' must be square (got %d x %d)", n, ncol(D)
        ))
    }

    K <- as.integer(K)
    if (length(K) != 1L || is.na(K) || K < 1L || K >= n) {
        stop(sprintf(
            "knn_from_matrix: K must be in [1, N-1]; got K=%d, N=%d", K, n
        ))
    }

    # ---- Default groups ----------------------------------------------------
    if (is.null(agroups)) agroups <- .default_groups(n)
    if (is.null(bgroups)) bgroups <- .default_groups(n)

    agroups <- as.integer(agroups)
    bgroups <- as.integer(bgroups)

    if (length(agroups) != n) {
        stop(sprintf(
            "knn_from_matrix: agroups length (%d) != nrow(D) (%d)",
            length(agroups), n
        ))
    }
    if (length(bgroups) != n) {
        stop(sprintf(
            "knn_from_matrix: bgroups length (%d) != nrow(D) (%d)",
            length(bgroups), n
        ))
    }

    rcpp_knn_from_distance_matrix(D, K, agroups, bgroups, isTRUE(sort_nbrs))
}


# ---------------------------------------------------------------------------
# knn_from_pca
# ---------------------------------------------------------------------------

#' K-nearest-neighbors from a PCA embedding matrix
#'
#' Computes K nearest neighbors by Euclidean distance from an N x D PCA
#' (or other embedding) matrix, without constructing the full N x N distance
#' matrix.  Optional same-group masking applies.
#'
#' @param pca_matrix Numeric matrix (N x D). Rows are samples, columns are
#'   embedding dimensions.
#' @param K Integer. Number of nearest neighbors. Must satisfy
#'   \code{1 <= K <= N - 1}.
#' @param agroups Integer vector of length N, or \code{NULL}. \code{NULL}
#'   assigns unique groups (no masking).
#' @param bgroups Integer vector of length N, or \code{NULL}.
#' @param sort_nbrs Logical. If \code{TRUE} (default), sort K neighbors by
#'   ascending Euclidean distance.
#' @return A \code{list} with two elements:
#'   \describe{
#'     \item{\code{knn_indices}}{Integer matrix (N x K). 1-based indices.}
#'     \item{\code{knn_distances}}{Numeric matrix (N x K). Euclidean distances.}
#'   }
#' @examples
#' \donttest{
#' pca <- matrix(rnorm(30), nrow = 10, ncol = 3)
#' knn <- knn_from_pca(pca, K = 3L)
#' }
#' @seealso \code{\link{knn_from_matrix}}, \code{\link{compute_tcrdist_kernel_pca}}
#' @export
knn_from_pca <- function(pca_matrix, K,
                          agroups   = NULL,
                          bgroups   = NULL,
                          sort_nbrs = TRUE) {

    # ---- Input validation ---------------------------------------------------
    if (!is.matrix(pca_matrix) || !is.numeric(pca_matrix)) {
        stop("knn_from_pca: 'pca_matrix' must be a numeric matrix")
    }
    n <- nrow(pca_matrix)

    K <- as.integer(K)
    if (length(K) != 1L || is.na(K) || K < 1L || K >= n) {
        stop(sprintf(
            "knn_from_pca: K must be in [1, N-1]; got K=%d, N=%d", K, n
        ))
    }

    # ---- Default groups ----------------------------------------------------
    if (is.null(agroups)) agroups <- .default_groups(n)
    if (is.null(bgroups)) bgroups <- .default_groups(n)

    agroups <- as.integer(agroups)
    bgroups <- as.integer(bgroups)

    if (length(agroups) != n) {
        stop(sprintf(
            "knn_from_pca: agroups length (%d) != nrow(pca_matrix) (%d)",
            length(agroups), n
        ))
    }
    if (length(bgroups) != n) {
        stop(sprintf(
            "knn_from_pca: bgroups length (%d) != nrow(pca_matrix) (%d)",
            length(bgroups), n
        ))
    }

    rcpp_knn_from_pca_matrix(pca_matrix, K, agroups, bgroups, isTRUE(sort_nbrs))
}
