# Kernel PCA on TCRdist distance matrices.
#
# Ported from Python CoNGA: conga/preprocess.py (lines 2203-2256).
# Matches scipy.linalg.eigh (via sklearn.decomposition.KernelPCA):
#   1. Build a Gram (kernel) matrix from the distance matrix.
#   2. Double-center the Gram matrix.
#   3. Eigendecompose and keep top positive-eigenvalue components.
#   4. Scale eigenvectors by sqrt(eigenvalues) to produce embeddings.
#
# R's base::eigen(symmetric=TRUE) uses the same LAPACK dsyevr routine as
# scipy.linalg.eigh, producing numerically identical results (~1e-10).


# ---------------------------------------------------------------------------
# .center_kernel_matrix (internal)
# ---------------------------------------------------------------------------

#' Double-center a kernel (Gram) matrix
#'
#' Applies the standard double-centering transform used by kernel PCA
#' (matches \code{sklearn.preprocessing.KernelCenterer}):
#' \deqn{K_c = K - 1_n K - K 1_n + 1_n K 1_n}
#' where \eqn{1_n} is an n-by-n matrix of \eqn{1/n}.
#'
#' @param K Numeric square matrix (the kernel / Gram matrix).
#' @return The double-centered kernel matrix (same dimensions as \code{K}).
#' @keywords internal
.center_kernel_matrix <- function(K) {
    n <- nrow(K)
    if (n == 0L) return(K)

    row_means <- rowMeans(K)
    col_means <- colMeans(K)
    grand_mean <- mean(K)

    K_c <- sweep(K, 1L, row_means, `-`)
    K_c <- sweep(K_c, 2L, col_means, `-`)
    K_c <- K_c + grand_mean

    K_c
}


# ---------------------------------------------------------------------------
# compute_tcrdist_kernel_pca
# ---------------------------------------------------------------------------

#' Kernel PCA on TCRdist distances
#'
#' Computes a low-dimensional embedding of TCRs by applying kernel PCA to
#' the pairwise TCRdist distance matrix.
#'
#' This implementation matches \code{scipy.linalg.eigh} (via sklearn's
#' \code{KernelPCA(kernel='precomputed')}). Both R's \code{base::eigen()}
#' and scipy use the same LAPACK \code{dsyevr} routine, so results are
#' numerically identical to ~1e-10 tolerance.
#'
#' Two kernel choices are supported:
#' \describe{
#'   \item{Default (\code{kernel = NULL})}{Linear kernel:
#'     \code{gram = pmax(0, 1 - D / Dmax)} where
#'     \code{Dmax = force_Dmax \%||\% max(D)}.}
#'   \item{Gaussian (\code{kernel = "gaussian"})}{RBF kernel:
#'     \code{gram = exp(-0.5 * (D / sdev)^2)}.}
#' }
#'
#' @param tcr_df A \code{data.frame} with at least columns \code{va},
#'   \code{cdr3a}, \code{vb}, \code{cdr3b}.
#' @param organism Character string. Organism key, e.g. \code{"human"} or
#'   \code{"mouse"}.
#' @param n_components Integer. Maximum number of PCA components to return.
#'   Clamped to \code{nrow(tcr_df)}. Default \code{50L}.
#' @param kernel \code{NULL} (default linear kernel) or \code{"gaussian"}.
#' @param gaussian_kernel_sdev Numeric. Standard deviation parameter for the
#'   Gaussian kernel. Ignored unless \code{kernel = "gaussian"}.
#'   Default \code{100}.
#' @param force_Dmax Numeric or \code{NULL}. If non-\code{NULL}, use this
#'   value instead of \code{max(D)} when computing the default kernel.
#'   Ignored when \code{kernel = "gaussian"}.
#' @param method Character. Eigen-decomposition method: \code{"auto"}
#'   (default, uses \code{RSpectra::eigs_sym()} when available for partial
#'   decomposition, falling back to \code{base::eigen()}), \code{"eigen"}
#'   (always uses \code{base::eigen()}, same LAPACK as scipy.linalg.eigh),
#'   or \code{"RSpectra"} (always uses \code{RSpectra::eigs_sym()}, same
#'   ARPACK as scipy.sparse.linalg.eigsh).
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{embeddings}}{Numeric matrix of dimensions
#'       N x n_components.}
#'     \item{\code{eigenvalues}}{Numeric vector of retained positive
#'       eigenvalues (decreasing order).}
#'     \item{\code{n_components}}{Integer. Number of components actually
#'       returned.}
#'   }
#' @examples
#' \donttest{
#' tcrs <- data.frame(
#'     va    = c("TRAV1-1*01", "TRAV1-2*01", "TRAV1-1*01"),
#'     cdr3a = c("CAVRDSSYKLIF", "CAVRDSNYQLIW", "CAVRDSSYKLIF"),
#'     vb    = c("TRBV19*01", "TRBV28*01", "TRBV19*01"),
#'     cdr3b = c("CASSIRSSYEQYF", "CASSLGQAYEQYF", "CASSIRSYEQYF"),
#'     stringsAsFactors = FALSE
#' )
#' result <- compute_tcrdist_kernel_pca(tcrs, "human", n_components = 2L)
#' str(result)
#' }
#' @seealso \code{\link{plot_tcr_scatter}}, \code{\link{knn_from_pca}}, \code{\link{tcrdist_matrix}}
#' @export
compute_tcrdist_kernel_pca <- function(tcr_df,
                                       organism,
                                       n_components = 50L,
                                       kernel = NULL,
                                       gaussian_kernel_sdev = 100,
                                       force_Dmax = NULL,
                                       method = c("auto", "eigen",
                                                   "RSpectra")) {
    # ---- Input validation --------------------------------------------------
    if (!is.data.frame(tcr_df)) {
        tcr_df <- as.data.frame(tcr_df, stringsAsFactors = FALSE)
    }
    required_cols <- c("va", "cdr3a", "vb", "cdr3b")
    missing_cols <- setdiff(required_cols, colnames(tcr_df))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "compute_tcrdist_kernel_pca: missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }
    if (!is.numeric(n_components) || length(n_components) != 1L ||
        n_components < 1L) {
        stop("compute_tcrdist_kernel_pca: 'n_components' must be a ",
             "positive integer")
    }
    n_components <- as.integer(n_components)
    if (!is.null(kernel) && !identical(kernel, "gaussian")) {
        stop(sprintf(
            "compute_tcrdist_kernel_pca: unrecognized kernel '%s'; use NULL or 'gaussian'",
            kernel
        ))
    }

    method <- match.arg(method)
    n <- nrow(tcr_df)
    if (n == 0L) {
        return(list(
            embeddings   = matrix(numeric(0L), nrow = 0L, ncol = 0L),
            eigenvalues  = numeric(0L),
            n_components = 0L
        ))
    }

    # "auto" selects RSpectra for partial decomposition when available
    if (method == "auto") {
        method <- if (requireNamespace("RSpectra", quietly = TRUE) &&
                      n_components < n) "RSpectra" else "eigen"
    }

    n_components <- min(n_components, n)

    # ---- Step 1: Compute TCRdist distance matrix ---------------------------
    if (n > 20000L) {
        warning(sprintf(
            "compute_tcrdist_kernel_pca: N=%d requires a full %d x %d distance matrix (%.1f GB). Consider reducing the dataset.",
            n, n, n, as.double(n) * n * 8 / 1e9
        ))
    }
    D <- tcrdist_matrix(tcr_df, organism)

    # ---- Step 2: Build Gram (kernel) matrix --------------------------------
    if (is.null(kernel)) {
        if (is.null(force_Dmax)) {
            force_Dmax <- max(D)
        }
        gram <- 1 - D / force_Dmax
        gram[gram < 0] <- 0
    } else {
        gram <- exp(-0.5 * (D / gaussian_kernel_sdev)^2)
    }

    # ---- Step 3: Double-center the Gram matrix -----------------------------
    gram_centered <- .center_kernel_matrix(gram)

    # Force symmetry (numerical safety for eigen())
    gram_centered <- (gram_centered + t(gram_centered)) / 2

    # ---- Step 4: Eigendecompose --------------------------------------------
    if (method == "eigen") {
        eig <- eigen(gram_centered, symmetric = TRUE)
        all_values  <- eig$values
        all_vectors <- eig$vectors
    } else {
        # RSpectra
        if (!requireNamespace("RSpectra", quietly = TRUE)) {
            stop("Package 'RSpectra' is required for method='RSpectra' ",
                 "but not installed.")
        }
        k_request <- min(n_components + 10L, n - 1L)
        eig <- RSpectra::eigs_sym(gram_centered, k = k_request,
                                  which = "LM")
        all_values  <- eig$values
        all_vectors <- eig$vectors
    }

    # ---- Step 5: Keep only positive eigenvalues ----------------------------
    pos_idx <- which(all_values > 0)
    if (length(pos_idx) == 0L) {
        return(list(
            embeddings   = matrix(0, nrow = n, ncol = 0L),
            eigenvalues  = numeric(0L),
            n_components = 0L
        ))
    }

    keep <- pos_idx[seq_len(min(n_components, length(pos_idx)))]
    eigenvalues <- all_values[keep]
    eigenvectors <- all_vectors[, keep, drop = FALSE]

    # ---- Step 6: Scale embeddings ------------------------------------------
    embeddings <- sweep(eigenvectors, 2L, sqrt(eigenvalues), `*`)

    k <- length(eigenvalues)

    list(
        embeddings   = embeddings,
        eigenvalues  = eigenvalues,
        n_components = k
    )
}
