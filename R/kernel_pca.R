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
#' @details For a symmetric \code{K} (the normal case: the Gram matrix is
#'   symmetrised before centering) \code{rowMeans == colMeans}, so the
#'   transform collapses to a single fused expression that avoids the two
#'   \code{sweep()} copies. A general (row/col separate) fallback is used when
#'   \code{K} is not symmetric.
#'
#' @param K Numeric square matrix (the kernel / Gram matrix).
#' @return The double-centered kernel matrix (same dimensions as \code{K}).
#' @keywords internal
.center_kernel_matrix <- function(K) {
    n <- nrow(K)
    if (n == 0L) return(K)

    if (isSymmetric(K)) {
        row_means <- rowMeans(K)
        return(K - outer(row_means, row_means, `+`) + mean(row_means))
    }

    row_means <- rowMeans(K)
    col_means <- colMeans(K)
    grand_mean <- mean(K)

    K_c <- sweep(K, 1L, row_means, `-`)
    K_c <- sweep(K_c, 2L, col_means, `-`)
    K_c <- K_c + grand_mean

    K_c
}


# ---------------------------------------------------------------------------
# RSpectra helpers (internal)
# ---------------------------------------------------------------------------

#' Is RSpectra usable in this session?
#' @keywords internal
.rspectra_available <- function() {
    requireNamespace("RSpectra", quietly = TRUE)
}

#' Matrix-free double-centered-kernel operator for ARPACK
#'
#' Computes \code{Kc \%*\% x} for the double-centered kernel \code{Kc} WITHOUT
#' materialising \code{Kc}, using (for symmetric \code{K}):
#' \deqn{Kc v = Kv - mean(Kv) - (rowSums(K)/n)(1'v) + (sum(K)/n^2)(1'v).}
#' \code{args} carries \code{K}, precomputed \code{row_sums} and \code{total}
#' (\code{= sum(K)}), and \code{n}, so per-iteration cost is one BLAS
#' \code{dgemv} plus O(n) work.
#' @keywords internal
.centered_kernel_matvec <- function(x, args) {
    Kv <- as.vector(args$K %*% x)
    sx <- sum(x)
    Kv - mean(Kv) - (args$row_sums / args$n) * sx + (args$total / args$n^2) * sx
}

#' Top-k eigenpairs of the double-centered kernel via ARPACK
#'
#' @param gram Symmetric raw kernel (Gram) matrix (not yet centered).
#' @param k_request Number of eigenpairs to request.
#' @param matrix_free If \code{TRUE}, apply centering implicitly via
#'   \code{.centered_kernel_matvec} (no centered matrix materialised); else
#'   double-centre \code{gram} densely first. Both use \code{which = "LA"}
#'   (largest algebraic) because the default kernel yields an indefinite
#'   centered Gram and only positive eigenvalues are wanted.
#' @return The \code{RSpectra::eigs_sym} result list (\code{values},
#'   \code{vectors}, \code{nconv}). Retries once with a larger Krylov subspace
#'   if the first solve under-converges.
#' @keywords internal
.kpca_eigs_rspectra <- function(gram, k_request, matrix_free = TRUE) {
    n <- nrow(gram)
    solve_once <- function(opts) {
        if (matrix_free) {
            row_sums <- rowSums(gram)
            args <- list(K = gram, row_sums = row_sums,
                         total = sum(row_sums), n = n)
            eigs_sym(.centered_kernel_matvec, k = k_request,
                     which = "LA", n = n, args = args, opts = opts)
        } else {
            eigs_sym(.center_kernel_matrix(gram), k = k_request,
                     which = "LA", opts = opts)
        }
    }
    eig <- solve_once(list(tol = 1e-12))
    nconv <- if (is.null(eig$nconv)) 0L else eig$nconv
    if (nconv < k_request) {
        eig <- solve_once(list(tol = 1e-12,
                               ncv = min(n, 4L * k_request + 1L),
                               maxitr = 3000L))
    }
    eig
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
#' @details
#' \strong{Performance.} With the default \code{method = "auto"}, when
#' \code{n_components} is small relative to \code{n} the decomposition uses a
#' matrix-free ARPACK partial solve (\code{RSpectra::eigs_sym}) that computes
#' only the requested components and never materialises the centered Gram
#' matrix -- typically an order of magnitude faster than a full
#' \code{base::eigen()} for large repertoires. Set \code{method = "eigen"} for a
#' bit-exact full LAPACK decomposition (e.g. to reproduce
#' \code{scipy.linalg.eigh}).
#'
#' The full-eigen path is bound by R's BLAS/LAPACK. On the reference
#' (unoptimised) BLAS an \code{n = 5000} decomposition can take minutes; linking
#' R against an optimised BLAS speeds all matrix math several-fold. On macOS this
#' is a one-line symlink to the Accelerate shim that already ships with R
#' (\code{libRblas.dylib} -> \code{libRblas.vecLib.dylib}); on Linux use
#' \code{update-alternatives} to select OpenBLAS. (This swaps BLAS only, not
#' LAPACK, but \code{dsyevr}'s dominant cost is BLAS-3, so most of the win
#' carries over.) The \code{"auto"}/ARPACK path sidesteps this by doing far less
#' arithmetic in the first place.
#'
#' @param tcr_df A \code{data.frame} with at least columns \code{va},
#'   \code{cdr3a}, \code{vb}, \code{cdr3b}. Optional if \code{dist_matrix}
#'   is provided.
#' @param organism Character string. Organism key, e.g. \code{"human"} or
#'   \code{"mouse"}. Optional if \code{dist_matrix} is provided.
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
#' @param dist_matrix Optional precomputed distance matrix. If provided,
#'   \code{tcr_df} and \code{organism} are not used for distance computation.
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
#' @importFrom RSpectra eigs_sym
#' @export
compute_tcrdist_kernel_pca <- function(tcr_df = NULL,
                                       organism = NULL,
                                       n_components = 50L,
                                       kernel = NULL,
                                       gaussian_kernel_sdev = 100,
                                       force_Dmax = NULL,
                                       method = c("auto", "eigen",
                                                   "RSpectra"),
                                       dist_matrix = NULL) {
    # ---- Input validation --------------------------------------------------
    if (!is.null(tcr_df) && !is.data.frame(tcr_df)) {
        tcr_df <- as.data.frame(tcr_df, stringsAsFactors = FALSE)
    }
    if (is.null(dist_matrix) && !is.null(tcr_df)) {
        required_cols <- c("va", "cdr3a", "vb", "cdr3b")
        missing_cols <- setdiff(required_cols, colnames(tcr_df))
        if (length(missing_cols) > 0L) {
            stop(sprintf(
                "compute_tcrdist_kernel_pca: missing required columns: %s",
                paste(missing_cols, collapse = ", ")
            ))
        }
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
    n <- if (!is.null(dist_matrix)) nrow(dist_matrix) else nrow(tcr_df)
    if (is.null(n) || n == 0L) {
        return(list(
            embeddings   = matrix(numeric(0L), nrow = 0L, ncol = 0L),
            eigenvalues  = numeric(0L),
            n_components = 0L
        ))
    }

    n_components <- min(n_components, n)

    # "auto" routes to the RSpectra (ARPACK) partial solver when it is a clear
    # win: enough clones that base::eigen()'s full O(n^3) solve dominates, and
    # few enough components (k <= n/2) that partial decomposition pays off.
    # Otherwise a single dense LAPACK call is faster and more robust.
    if (method == "auto") {
        use_partial <- .rspectra_available() && n >= 100L &&
                       n_components <= n %/% 2L && n_components <= n - 2L
        method <- if (use_partial) "RSpectra" else "eigen"
    }

    # ---- Step 1: Compute TCRdist distance matrix ---------------------------
    if (is.null(dist_matrix) && n > 20000L) {
        warning(sprintf(
            "compute_tcrdist_kernel_pca: N=%d requires a full %d x %d distance matrix (%.1f GB). Consider reducing the dataset.",
            n, n, n, as.double(n) * n * 8 / 1e9
        ))
    }
    D <- .get_dist_matrix(tcr_df, organism, dist_matrix)

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
    # D is no longer needed; free it before the eigensolve when we own it.
    if (is.null(dist_matrix)) rm(D)

    # ---- Step 3: Symmetrise (guard for user-supplied dist_matrix) ----------
    # Both the fused centering and the matrix-free operator assume a symmetric
    # Gram matrix. tcrdist_matrix() is symmetric; a user-supplied dist_matrix
    # may not be. Symmetrising here (before centering, which preserves it)
    # replaces the previous post-centering re-symmetrisation.
    if (!isSymmetric(gram)) {
        gram <- (gram + t(gram)) / 2
    }

    # ---- Step 4: Eigendecompose --------------------------------------------
    if (method == "RSpectra") {
        if (!.rspectra_available()) {
            stop("Package 'RSpectra' is required for method='RSpectra' ",
                 "but not installed.")
        }
        k_request <- min(n_components, n - 1L)
        eig <- .kpca_eigs_rspectra(gram, k_request, matrix_free = TRUE)
        nconv <- if (is.null(eig$nconv)) 0L else eig$nconv
        if (nconv < k_request) {
            warning(sprintf(
                "compute_tcrdist_kernel_pca: ARPACK converged %d/%d eigenpairs; falling back to base::eigen().",
                nconv, k_request))
            method <- "eigen"
        } else {
            all_values  <- eig$values
            all_vectors <- eig$vectors
        }
    }
    if (method == "eigen") {
        eig <- eigen(.center_kernel_matrix(gram), symmetric = TRUE)
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
