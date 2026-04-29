#' Compute UMAP embedding from TCRdist kernel PCA
#'
#' Reduces kernel PCA embeddings (or raw TCR data) to a low-dimensional UMAP
#' representation suitable for visualization with \code{\link{plot_tcr_scatter}}.
#' Uses the \pkg{uwot} package for the UMAP computation.
#'
#' Two input modes are supported:
#' \enumerate{
#'   \item \strong{From pre-computed PCA}: supply \code{pca_embeddings} (an
#'     N x D matrix, e.g. from \code{compute_tcrdist_kernel_pca()$embeddings}).
#'   \item \strong{From raw TCR data}: supply \code{tcr_df} and
#'     \code{organism}. Kernel PCA is computed internally with
#'     \code{n_components_pca} components before UMAP.
#' }
#'
#' @param tcr_df Data.frame with TCR columns (\code{va}, \code{cdr3a},
#'   \code{vb}, \code{cdr3b}). Used only if \code{pca_embeddings} is
#'   \code{NULL}.
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#'   Required when \code{tcr_df} is used.
#' @param pca_embeddings Numeric matrix (N x D). Pre-computed kernel PCA
#'   embeddings. If provided, \code{tcr_df} and \code{organism} are ignored.
#' @param n_components_pca Integer. Number of kernel PCA components to compute
#'   when using the \code{tcr_df} input path. Default \code{50L}.
#' @param n_components Integer. Number of UMAP output dimensions. Default
#'   \code{2L}.
#' @param n_neighbors Integer. Size of local neighborhood for UMAP manifold
#'   approximation. Default \code{15L}.
#' @param min_dist Numeric. Minimum distance between embedded points. Controls
#'   how tightly UMAP packs points together. Default \code{0.1}.
#' @param metric Character string. Distance metric for UMAP neighbor search.
#'   Default \code{"euclidean"}.
#' @param n_threads Integer. Number of threads for neighbor search and
#'   optimization. Default \code{1L}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility. If
#'   non-\code{NULL}, \code{set.seed()} is called before UMAP computation.
#'   Default \code{NULL}.
#' @param ... Additional arguments passed to \code{uwot::umap()}.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{embeddings}}{Numeric matrix (N x \code{n_components}).
#'       UMAP coordinates.}
#'     \item{\code{pca_embeddings}}{Numeric matrix. The PCA input used.}
#'     \item{\code{n_components}}{Integer. Number of UMAP dimensions.}
#'   }
#'
#' @examples
#' \donttest{
#' data(dash)
#' # From raw data (computes kernel PCA internally)
#' umap <- compute_tcrdist_umap(dash[1:100, ], "mouse", seed = 42)
#' dim(umap$embeddings)  # 100 x 2
#'
#' # From pre-computed PCA
#' pca <- compute_tcrdist_kernel_pca(dash[1:100, ], "mouse", n_components = 20L)
#' umap <- compute_tcrdist_umap(pca_embeddings = pca$embeddings, seed = 42)
#' }
#'
#' @seealso \code{\link{compute_tcrdist_kernel_pca}},
#'   \code{\link{plot_tcr_scatter}}
#' @export
compute_tcrdist_umap <- function(tcr_df = NULL,
                                  organism = NULL,
                                  pca_embeddings = NULL,
                                  n_components_pca = 50L,
                                  n_components = 2L,
                                  n_neighbors = 15L,
                                  min_dist = 0.1,
                                  metric = "euclidean",
                                  n_threads = 1L,
                                  seed = NULL,
                                  ...) {
    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required for compute_tcrdist_umap(). ",
             "Install it with: install.packages(\"uwot\")",
             call. = FALSE)
    }

    # ---- Resolve PCA embeddings -----------------------------------------------
    if (is.null(pca_embeddings)) {
        if (is.null(tcr_df)) {
            stop("compute_tcrdist_umap: provide either 'pca_embeddings' or ",
                 "'tcr_df' + 'organism'", call. = FALSE)
        }
        if (is.null(organism)) {
            stop("compute_tcrdist_umap: 'organism' is required when using ",
                 "'tcr_df' input", call. = FALSE)
        }
        pca_result <- compute_tcrdist_kernel_pca(
            tcr_df, organism,
            n_components = as.integer(n_components_pca)
        )
        pca_embeddings <- pca_result$embeddings
    }

    pca_embeddings <- as.matrix(pca_embeddings)
    n <- nrow(pca_embeddings)
    if (n < 2L) {
        stop("compute_tcrdist_umap: need at least 2 samples", call. = FALSE)
    }

    # Clamp n_neighbors to dataset size
    n_neighbors <- min(as.integer(n_neighbors), n - 1L)

    # ---- Run UMAP -------------------------------------------------------------
    if (!is.null(seed)) set.seed(seed)

    umap_coords <- uwot::umap(
        pca_embeddings,
        n_neighbors  = n_neighbors,
        n_components = as.integer(n_components),
        min_dist     = min_dist,
        metric       = metric,
        n_threads    = as.integer(n_threads),
        ...
    )

    list(
        embeddings     = umap_coords,
        pca_embeddings = pca_embeddings,
        n_components   = as.integer(n_components)
    )
}
