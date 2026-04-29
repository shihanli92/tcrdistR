# ---------------------------------------------------------------------------
# Internal helpers ported from rconga/R/clustering.R
# ---------------------------------------------------------------------------

#' Compute per-point bandwidth for fuzzy simplicial set
#'
#' Binary search for sigma (per-point bandwidth) and rho (local connectivity
#' distance) used by UMAP's fuzzy simplicial set construction.  Port of
#' rconga's \code{.smooth_knn_dist()}.
#'
#' @param knn_distances Numeric matrix (N x K).  Sorted KNN distances.
#' @param local_connectivity Numeric scalar.  Default \code{1.0}.
#' @param bandwidth Numeric scalar.  Default \code{1.0}.
#' @param n_iter Integer.  Max binary search iterations.  Default \code{64L}.
#' @return Named list with \code{rho} and \code{sigma} (numeric vectors of
#'   length N).
#' @keywords internal
.smooth_knn_dist <- function(knn_distances,
                              local_connectivity = 1.0,
                              bandwidth = 1.0,
                              n_iter = 64L) {
    n <- nrow(knn_distances)
    k <- ncol(knn_distances)
    target <- log2(k + 1L) * bandwidth

    rho   <- numeric(n)
    sigma <- numeric(n)

    SMOOTH_K_TOLERANCE <- 1e-5
    MIN_K_DIST_SCALE   <- 1e-3

    # ---- Compute rho (distance to nearest non-identical neighbor) ------------
    index <- floor(local_connectivity)
    interpolation <- local_connectivity - index

    for (i in seq_len(n)) {
        non_zero <- knn_distances[i, knn_distances[i, ] > 0]
        if (length(non_zero) >= local_connectivity) {
            if (index > 0L) {
                rho[i] <- non_zero[index]
                if (interpolation > SMOOTH_K_TOLERANCE) {
                    rho[i] <- rho[i] +
                        interpolation * (non_zero[index + 1L] - non_zero[index])
                }
            } else {
                rho[i] <- interpolation * non_zero[1L]
            }
        } else if (length(non_zero) > 0L) {
            rho[i] <- max(non_zero)
        }
    }

    # ---- Binary search for sigma ---------------------------------------------
    mean_all_dist <- mean(knn_distances)

    for (i in seq_len(n)) {
        lo <- 0.0
        hi <- Inf
        mid <- 1.0
        dists_i <- knn_distances[i, ]

        for (iter in seq_len(n_iter)) {
            d_shifted <- dists_i - rho[i]
            psum <- sum(ifelse(d_shifted > 0, exp(-d_shifted / mid), 1.0))

            if (abs(psum - target) < SMOOTH_K_TOLERANCE) break

            if (psum > target) {
                hi  <- mid
                mid <- (lo + hi) / 2.0
            } else {
                lo <- mid
                if (is.infinite(hi)) {
                    mid <- mid * 2.0
                } else {
                    mid <- (lo + hi) / 2.0
                }
            }
        }

        sigma[i] <- mid

        mean_d <- if (rho[i] > 0) mean(dists_i) else mean_all_dist
        if (sigma[i] < MIN_K_DIST_SCALE * mean_d) {
            sigma[i] <- MIN_K_DIST_SCALE * mean_d
        }
    }

    list(rho = rho, sigma = sigma)
}


# ---------------------------------------------------------------------------
# .build_knn_graph
# ---------------------------------------------------------------------------

#' Build a fuzzy simplicial set graph from KNN data
#'
#' Constructs a sparse N x N adjacency matrix using UMAP's fuzzy simplicial
#' set algorithm.  Edge weights are \code{exp(-(d - rho) / sigma)} and the
#' graph is symmetrized via fuzzy union: \code{W + W^T - W * W^T}.  Port of
#' rconga's \code{.build_knn_graph()}, adapted for 1-based indices.
#'
#' @param knn_indices Integer matrix (N x K).  1-based neighbor indices (as
#'   returned by \code{\link{tcrdist_knn}}).
#' @param knn_distances Numeric matrix (N x K).  Corresponding distances.
#' @param n_cells Integer scalar.  Total number of observations.
#' @return Sparse \code{dgCMatrix} (N x N) with symmetrized fuzzy membership
#'   weights.
#' @keywords internal
.build_knn_graph <- function(knn_indices, knn_distances, n_cells) {
    n <- nrow(knn_indices)
    k <- ncol(knn_indices)

    smooth <- .smooth_knn_dist(knn_distances)

    # knn_indices are 1-based (tcrdistR convention)
    rows <- as.vector(row(knn_indices))
    cols <- as.vector(knn_indices)

    dist_flat      <- as.vector(knn_distances)
    rho_expanded   <- rep(smooth$rho, times = k)
    sigma_expanded <- rep(smooth$sigma, times = k)

    d_shifted <- dist_flat - rho_expanded
    weights <- ifelse(
        d_shifted <= 0 | sigma_expanded == 0,
        1.0,
        exp(-d_shifted / sigma_expanded)
    )

    # Zero out self-loops
    self_loop <- (rows == cols)
    weights[self_loop] <- 0.0

    W <- Matrix::sparseMatrix(
        i    = rows,
        j    = cols,
        x    = weights,
        dims = c(n_cells, n_cells)
    )

    # Symmetrize: fuzzy set union W + W^T - W * W^T
    Wt <- Matrix::t(W)
    W + Wt - W * Wt
}


# ---------------------------------------------------------------------------
# .run_umap_from_knn
# ---------------------------------------------------------------------------

#' Run UMAP from precomputed K-nearest-neighbor data
#'
#' Passes precomputed neighbor indices and distances to
#' \code{\link[uwot]{umap}} via its \code{nn_method} list interface.  Port of
#' rconga's \code{.run_umap_from_knn()}, adapted for 1-based indices.
#'
#' @param knn_indices Integer matrix (N x K).  1-based neighbor indices.
#' @param knn_distances Numeric matrix (N x K).  Distances.
#' @param n_components Integer.  UMAP output dimensions.  Default \code{2L}.
#' @param min_dist Numeric.  UMAP min_dist.  Default \code{0.5}.
#' @param spread Numeric.  UMAP spread.  Default \code{1.0}.
#' @param seed Integer or \code{NULL}.  Random seed.
#' @param n_threads Integer.  Threads for optimization.  Default \code{1L}.
#' @return Numeric matrix (N x n_components).
#' @keywords internal
.run_umap_from_knn <- function(knn_indices,
                                knn_distances,
                                n_components = 2L,
                                min_dist = 0.5,
                                spread = 1.0,
                                seed = NULL,
                                n_threads = 1L) {
    # uwot expects 1-based indices --- tcrdistR already uses 1-based
    idx <- knn_indices
    storage.mode(idx) <- "integer"

    if (!is.null(seed)) set.seed(seed)

    uwot::umap(
        X            = knn_distances,
        nn_method    = list(idx = idx, dist = knn_distances),
        n_components = as.integer(n_components),
        min_dist     = min_dist,
        spread       = spread,
        n_threads    = as.integer(n_threads),
        ret_nn       = FALSE
    )
}


# ---------------------------------------------------------------------------
# .run_clustering
# ---------------------------------------------------------------------------

#' Run graph-based clustering on a fuzzy simplicial set graph
#'
#' Performs Leiden or Louvain community detection on the KNN graph.  When
#' \code{method} is \code{NULL}, Leiden is tried first with Louvain as
#' fallback.  Port of rconga's \code{.run_clustering()}.
#'
#' @param knn_graph Sparse matrix (N x N) from \code{.build_knn_graph()}.
#' @param resolution Numeric.  Resolution parameter.  Default \code{1.0}.
#' @param method Character or \code{NULL}.  \code{"leiden"}, \code{"louvain"},
#'   or \code{NULL} (try Leiden first).
#' @return Integer vector of length N with 0-based cluster IDs.
#' @keywords internal
.run_clustering <- function(knn_graph,
                             resolution = 1.0,
                             method = NULL) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required for clustering. ",
             "Install it with: install.packages(\"igraph\")",
             call. = FALSE)
    }

    graph <- igraph::graph_from_adjacency_matrix(
        knn_graph, mode = "undirected", weighted = TRUE, diag = FALSE
    )

    if (is.null(method)) {
        result <- tryCatch({
            igraph::cluster_leiden(graph, resolution = resolution,
                                   objective_function = "modularity")
        }, error = function(e) {
            message("Leiden clustering failed, falling back to Louvain: ",
                    conditionMessage(e))
            igraph::cluster_louvain(graph, resolution = resolution)
        })
    } else {
        method <- match.arg(method, c("leiden", "louvain"))
        if (method == "leiden") {
            result <- igraph::cluster_leiden(graph, resolution = resolution,
                                             objective_function = "modularity")
        } else {
            result <- igraph::cluster_louvain(graph, resolution = resolution)
        }
    }

    as.integer(igraph::membership(result)) - 1L
}


# ---------------------------------------------------------------------------
# .calc_nndists
# ---------------------------------------------------------------------------

#' Weighted nearest-neighbor distances
#'
#' Computes a weighted average of each row of the sorted KNN distance matrix.
#' Weights decrease linearly from 1.0 (nearest) to 1/K (farthest).  Port of
#' rconga's \code{.calc_nndists()}.
#'
#' @param knn_distances Numeric matrix (N x K).  Sorted ascending.
#' @return Numeric vector of length N.
#' @keywords internal
.calc_nndists <- function(knn_distances) {
    k <- ncol(knn_distances)
    if (nrow(knn_distances) == 0L || k == 0L) return(numeric(0L))

    wts <- seq(1.0, 1.0 / k, length.out = k)
    wts <- wts / sum(wts)
    as.numeric(knn_distances %*% wts)
}


# ---------------------------------------------------------------------------
# .umap_from_pca (simple Euclidean UMAP on PCA embeddings)
# ---------------------------------------------------------------------------

#' @keywords internal
.umap_from_pca <- function(pca_embeddings, n_components, n_neighbors,
                            min_dist, metric, n_threads, seed, ...) {
    pca_embeddings <- as.matrix(pca_embeddings)
    n <- nrow(pca_embeddings)
    if (n < 2L) {
        stop("compute_tcrdist_umap: need at least 2 samples", call. = FALSE)
    }
    n_neighbors <- min(as.integer(n_neighbors), n - 1L)

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
        n_components   = as.integer(n_components),
        method         = "pca"
    )
}


# ===========================================================================
# Public API
# ===========================================================================

#' Compute UMAP embedding from TCRdist data
#'
#' Two input modes are supported:
#'
#' \enumerate{
#'   \item \strong{KNN path} (recommended): supply \code{tcr_df} and
#'     \code{organism}.
#'     Computes TCRdist K-nearest-neighbors with group masking (clones sharing
#'     an identical alpha or beta chain are excluded from each other's
#'     neighborhoods), builds a fuzzy simplicial set graph, and runs UMAP from
#'     the precomputed KNN via \code{\link[uwot]{umap}}.  This preserves the
#'     TCRdist metric faithfully.
#'   \item \strong{PCA path}: supply \code{pca_embeddings} (an N x D matrix,
#'     e.g.\sspace{}from \code{compute_tcrdist_kernel_pca()$embeddings}).
#'     Runs standard UMAP in Euclidean space on the PCA coordinates.
#' }
#'
#' When \code{cluster = TRUE}, Leiden (or Louvain) community detection is
#' performed on the fuzzy KNN graph (KNN path only; requires \pkg{igraph}).
#'
#' @param tcr_df Data.frame with TCR columns (\code{va}, \code{cdr3a},
#'   \code{vb}, \code{cdr3b}; plus \code{ja}, \code{jb} for group masking).
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#'   Required when \code{tcr_df} is used.
#' @param pca_embeddings Numeric matrix (N x D).  Pre-computed kernel PCA
#'   embeddings.  If provided, the PCA path is used.
#' @param n_components Integer.  Number of UMAP output dimensions.  Default
#'   \code{2L}.
#' @param n_neighbors Integer.  Number of nearest neighbors.  Default
#'   \code{15L}.
#' @param min_dist Numeric.  UMAP min_dist.  Default \code{0.1}.
#' @param spread Numeric.  UMAP spread parameter.  Default \code{1.0}.
#' @param metric Character.  Distance metric for UMAP neighbor search (PCA
#'   path only).  Default \code{"euclidean"}.
#' @param n_threads Integer.  Threads for UMAP optimization.  Default
#'   \code{1L}.
#' @param seed Integer or \code{NULL}.  Random seed.  Default \code{NULL}.
#' @param cluster Logical.  If \code{TRUE}, run graph-based clustering on the
#'   KNN graph (KNN path only; requires \pkg{igraph}).  Default \code{FALSE}.
#' @param clustering_resolution Numeric.  Resolution for community detection.
#'   Default \code{1.0}.
#' @param clustering_method Character or \code{NULL}.  \code{"leiden"},
#'   \code{"louvain"}, or \code{NULL} (try Leiden first).  Default
#'   \code{NULL}.
#' @param ... Additional arguments passed to \code{\link[uwot]{umap}} (PCA
#'   path only).
#'
#' @return A named list.  Elements depend on the input path:
#'
#'   \strong{KNN path} (\code{tcr_df} + \code{organism}):
#'   \describe{
#'     \item{\code{embeddings}}{Numeric matrix (N x \code{n_components}).}
#'     \item{\code{knn_indices}}{Integer matrix (N x K).  1-based.}
#'     \item{\code{knn_distances}}{Numeric matrix (N x K).}
#'     \item{\code{knn_graph}}{Sparse \code{dgCMatrix} (N x N).  Fuzzy
#'       simplicial set.}
#'     \item{\code{clusters}}{Integer vector (0-based) or \code{NULL}.}
#'     \item{\code{nndists}}{Numeric vector.  Weighted NN distances.}
#'     \item{\code{n_neighbors}}{Final K (may have been increased for
#'       connectivity).}
#'     \item{\code{n_components}}{Integer.}
#'     \item{\code{method}}{\code{"knn"}.}
#'   }
#'
#'   \strong{PCA path} (\code{pca_embeddings}):
#'   \describe{
#'     \item{\code{embeddings}}{Numeric matrix (N x \code{n_components}).}
#'     \item{\code{pca_embeddings}}{The PCA input matrix.}
#'     \item{\code{n_components}}{Integer.}
#'     \item{\code{method}}{\code{"pca"}.}
#'   }
#'
#' @examples
#' \donttest{
#' data(dash)
#' sub <- dash[1:200, ]
#'
#' # KNN path (recommended): uses TCRdist directly
#' umap <- compute_tcrdist_umap(sub, "mouse", seed = 42)
#' dim(umap$embeddings)  # 200 x 2
#'
#' # With clustering
#' umap <- compute_tcrdist_umap(sub, "mouse", seed = 42, cluster = TRUE)
#' table(umap$clusters)
#'
#' # PCA path: from pre-computed kernel PCA
#' pca <- compute_tcrdist_kernel_pca(sub, "mouse", n_components = 20L)
#' umap <- compute_tcrdist_umap(pca_embeddings = pca$embeddings, seed = 42)
#' }
#'
#' @seealso \code{\link{tcrdist_knn}}, \code{\link{setup_tcr_groups}},
#'   \code{\link{compute_tcrdist_kernel_pca}}, \code{\link{plot_tcr_scatter}}
#' @export
compute_tcrdist_umap <- function(tcr_df = NULL,
                                  organism = NULL,
                                  pca_embeddings = NULL,
                                  n_components = 2L,
                                  n_neighbors = 15L,
                                  min_dist = 0.1,
                                  spread = 1.0,
                                  metric = "euclidean",
                                  n_threads = 1L,
                                  seed = NULL,
                                  cluster = FALSE,
                                  clustering_resolution = 1.0,
                                  clustering_method = NULL,
                                  ...) {
    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required for compute_tcrdist_umap(). ",
             "Install it with: install.packages(\"uwot\")",
             call. = FALSE)
    }

    # ---- PCA path: simple Euclidean UMAP on embeddings -----------------------
    if (!is.null(pca_embeddings)) {
        return(.umap_from_pca(pca_embeddings, n_components, n_neighbors,
                               min_dist, metric, n_threads, seed, ...))
    }

    # ---- KNN path: TCRdist KNN -> fuzzy graph -> UMAP ------------------------
    if (is.null(tcr_df)) {
        stop("compute_tcrdist_umap: provide either 'pca_embeddings' or ",
             "'tcr_df' + 'organism'", call. = FALSE)
    }
    if (is.null(organism)) {
        stop("compute_tcrdist_umap: 'organism' is required when using ",
             "'tcr_df' input", call. = FALSE)
    }

    n <- nrow(tcr_df)
    if (n < 2L) {
        stop("compute_tcrdist_umap: need at least 2 samples", call. = FALSE)
    }
    n_neighbors <- as.integer(n_neighbors)
    n_components <- as.integer(n_components)

    # Step 1: TCR groups for chain masking
    groups <- tryCatch(
        setup_tcr_groups(tcr_df),
        error = function(e) NULL
    )
    agroups <- groups$agroups
    bgroups <- groups$bgroups

    # Step 2: KNN with group masking + connectivity check
    knn <- NULL
    knn_graph <- NULL

    repeat {
        knn <- tcrdist_knn(
            tcr_df, organism, K = n_neighbors,
            agroups = agroups, bgroups = bgroups, sort_nbrs = TRUE
        )

        knn_graph <- .build_knn_graph(
            knn$knn_indices, knn$knn_distances, n
        )

        # Check connected components if igraph is available
        if (requireNamespace("igraph", quietly = TRUE)) {
            g_tmp <- igraph::graph_from_adjacency_matrix(
                knn_graph, mode = "undirected", weighted = TRUE, diag = FALSE
            )
            n_comp <- igraph::components(g_tmp)$no
            if (n_comp > 2L * n_components) {
                message(
                    "KNN graph has ", n_comp, " connected components; ",
                    "increasing n_neighbors from ", n_neighbors,
                    " to ", 2L * n_neighbors
                )
                n_neighbors <- 2L * n_neighbors
                next
            }
        }
        break
    }

    # Step 3: UMAP from precomputed KNN
    umap_coords <- .run_umap_from_knn(
        knn$knn_indices, knn$knn_distances,
        n_components = n_components, min_dist = min_dist,
        spread = spread, seed = seed, n_threads = n_threads
    )

    # Step 4: Optional graph-based clustering
    clusters <- NULL
    if (isTRUE(cluster)) {
        clusters <- .run_clustering(
            knn_graph,
            resolution = clustering_resolution,
            method = clustering_method
        )
    }

    # Step 5: Weighted nearest-neighbor distances
    nndists <- .calc_nndists(knn$knn_distances)

    list(
        embeddings    = umap_coords,
        knn_indices   = knn$knn_indices,
        knn_distances = knn$knn_distances,
        knn_graph     = knn_graph,
        clusters      = clusters,
        nndists       = nndists,
        n_neighbors   = n_neighbors,
        n_components  = n_components,
        method        = "knn"
    )
}
