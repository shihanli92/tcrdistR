# Hierarchical clustering and neighborhood-based statistical testing
#
# Ported from Python tcrdist3 rep_diff.py and tree.py.


# ---------------------------------------------------------------------------
# .neighborhood_tally  (internal)
# ---------------------------------------------------------------------------

#' Tally category membership within each TCR's neighborhood
#'
#' For each TCR, counts how many neighbors belong to each category in
#' \code{membership}. Neighborhood is defined by a distance matrix and
#' radius threshold.
#'
#' @param dist_mat Numeric matrix. Pairwise distance matrix (N x N).
#' @param membership Factor or character vector of length N.
#' @param radius Numeric. Maximum distance to be considered a neighbor.
#' @return A matrix: rows = TCRs, cols = categories, values = counts.
#' @keywords internal
#' @noRd
.neighborhood_tally <- function(dist_mat, membership, radius) {
    n <- nrow(dist_mat)
    membership <- as.factor(membership)
    lvls <- levels(membership)
    n_lvls <- length(lvls)

    tally <- matrix(0L, nrow = n, ncol = n_lvls,
                    dimnames = list(NULL, lvls))

    for (i in seq_len(n)) {
        nbrs <- which(dist_mat[i, ] <= radius)
        for (j in nbrs) {
            cat_j <- as.character(membership[j])
            tally[i, cat_j] <- tally[i, cat_j] + 1L
        }
    }

    tally
}


# ---------------------------------------------------------------------------
# tcrdist_hclust  (exported)
# ---------------------------------------------------------------------------

#' Hierarchical clustering of TCRs by TCRdist
#'
#' Convenience wrapper that computes the pairwise TCRdist matrix and
#' performs hierarchical clustering via \code{stats::hclust()}.
#'
#' @param tcr_df Data.frame with TCR columns (optional if \code{dist_matrix}
#'   is provided).
#' @param organism Character string (\code{"human"} or \code{"mouse"})
#'   (optional if \code{dist_matrix} is provided).
#' @param method Clustering method for \code{stats::hclust()}. Default
#'   \code{"average"} (UPGMA).
#' @param max_tcrs Integer. Subsample if N exceeds this. Default \code{2000L}.
#' @param dist_matrix Optional precomputed distance matrix. If provided,
#'   \code{tcr_df} and \code{organism} are not used for distance computation.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{hclust}}{An \code{hclust} object.}
#'     \item{\code{dist_matrix}}{The pairwise distance matrix used.}
#'     \item{\code{indices}}{Integer vector of row indices used (after
#'       potential subsampling).}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- tcrdist_hclust(tcr_df, "human")
#' plot(result$hclust)
#' clusters <- cutree(result$hclust, k = 5)
#' }
#'
#' @seealso \code{\link{cluster_tcrs}}, \code{\link{neighborhood_test}}, \code{\link{plot_tcrdist_dendrogram}}
#' @export
tcrdist_hclust <- function(tcr_df = NULL, organism = NULL, method = "average",
                            max_tcrs = 2000L, dist_matrix = NULL) {
    n <- if (!is.null(dist_matrix)) nrow(dist_matrix) else nrow(tcr_df)
    stopifnot(n >= 2L)

    indices <- seq_len(n)
    if (n > max_tcrs) {
        indices <- sort(sample(n, max_tcrs))
        if (!is.null(tcr_df)) tcr_df <- tcr_df[indices, , drop = FALSE]
        if (!is.null(dist_matrix)) dist_matrix <- dist_matrix[indices, indices, drop = FALSE]
    }

    dist_mat <- .get_dist_matrix(tcr_df, organism, dist_matrix)
    hc <- stats::hclust(stats::as.dist(dist_mat), method = method)

    list(
        hclust = hc,
        dist_matrix = dist_mat,
        indices = indices
    )
}


# ---------------------------------------------------------------------------
# .get_dist_matrix  (internal)
# ---------------------------------------------------------------------------

#' Get or compute a dense distance matrix
#'
#' Returns \code{dist_matrix} if provided, otherwise computes one from
#' \code{tcr_df} and \code{organism}.
#'
#' @param tcr_df Data.frame with TCR columns, or \code{NULL}.
#' @param organism Character string, or \code{NULL}.
#' @param dist_matrix Precomputed distance matrix, or \code{NULL}.
#' @return A numeric matrix (N x N).
#' @keywords internal
#' @noRd
.get_dist_matrix <- function(tcr_df, organism, dist_matrix) {
    if (!is.null(dist_matrix)) {
        return(as.matrix(dist_matrix))
    }
    if (is.null(tcr_df) || is.null(organism)) {
        stop("Either 'dist_matrix' or both 'tcr_df' and 'organism' must be ",
             "provided", call. = FALSE)
    }
    tcrdist_matrix(tcr_df, organism)
}


# ---------------------------------------------------------------------------
# .dbscan_auto_eps  (internal)
# ---------------------------------------------------------------------------

#' Auto-detect DBSCAN eps using k-distance knee
#'
#' Computes the sorted k-th nearest neighbor distances and finds the point
#' of maximum curvature (the "knee") as the recommended eps.
#'
#' @param dist_mat Numeric matrix (N x N). Distance matrix.
#' @param min_pts Integer. The minPts parameter for DBSCAN.
#' @return Numeric scalar. The auto-detected eps value.
#' @keywords internal
#' @noRd
.dbscan_auto_eps <- function(dist_mat, min_pts) {
    n <- nrow(dist_mat)
    k <- min(min_pts, n - 1L)

    # For each point, get sorted distances to all others and take the k-th
    k_dists <- numeric(n)
    for (i in seq_len(n)) {
        sorted_d <- sort(dist_mat[i, -i])
        k_dists[i] <- sorted_d[k]
    }

    # Sort k-distances in ascending order
    k_dists <- sort(k_dists)

    # Find the knee: point of maximum second derivative (discrete curvature)
    if (length(k_dists) < 3L) return(stats::median(k_dists))

    # Second differences approximate curvature
    d2 <- diff(diff(k_dists))
    knee_idx <- which.max(d2) + 1L
    eps <- k_dists[knee_idx]

    message("Auto-detected DBSCAN eps: ", round(eps, 1))
    eps
}


# ---------------------------------------------------------------------------
# cluster_tcrs  (exported)
# ---------------------------------------------------------------------------

#' Cluster TCRs using various algorithms
#'
#' Assigns TCR clonotypes to clusters using one of five methods:
#' hierarchical clustering, Leiden or Louvain community detection,
#' DBSCAN density-based clustering, or k-medoids (PAM).
#'
#' Distances can be precomputed via \code{dist_matrix} or computed
#' internally from \code{tcr_df} and \code{organism}.
#'
#' @param tcr_df Data.frame with TCR columns. Required unless
#'   \code{dist_matrix} is provided.
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#'   Required unless \code{dist_matrix} is provided.
#' @param method Clustering method. One of \code{"hierarchical"} (default),
#'   \code{"leiden"}, \code{"louvain"}, \code{"dbscan"}, or
#'   \code{"kmedoids"}.
#' @param dist_matrix Optional precomputed N x N distance matrix. If
#'   provided, \code{tcr_df} and \code{organism} are only needed for
#'   graph-based methods (Leiden/Louvain) when KNN must be computed.
#' @param k Integer. Number of clusters for \code{"hierarchical"} and
#'   \code{"kmedoids"}. For hierarchical, exactly one of \code{k} or
#'   \code{h} must be specified. For kmedoids, required.
#' @param h Numeric. Height for dendrogram cutting (hierarchical only).
#' @param hclust_method Character. Agglomeration method for
#'   \code{stats::hclust()}. Default \code{"average"} (UPGMA).
#' @param resolution Numeric. Resolution parameter for Leiden/Louvain.
#'   Higher values yield more clusters. Default \code{1.0}.
#' @param n_neighbors Integer. Number of nearest neighbors for KNN graph
#'   construction (Leiden/Louvain). Default \code{10L}.
#' @param eps Numeric or \code{NULL}. DBSCAN neighborhood radius. If
#'   \code{NULL}, auto-detected from the k-distance knee.
#' @param min_pts Integer. DBSCAN minimum points for core point.
#'   Default \code{5L}.
#'
#' @return An integer vector of cluster assignments (length N).
#'   \itemize{
#'     \item For \code{hierarchical}, \code{leiden}, \code{louvain}, and
#'       \code{kmedoids}: 1-based cluster IDs.
#'     \item For \code{dbscan}: 0-based (0 = noise/unassigned), with
#'       cluster IDs starting at 1.
#'   }
#'
#' @examples
#' \dontrun{
#' data(dash)
#' sub <- dash[1:50, ]
#'
#' # Hierarchical (default)
#' clusters <- cluster_tcrs(sub, "mouse", k = 5)
#'
#' # K-medoids
#' km <- cluster_tcrs(sub, "mouse", method = "kmedoids", k = 4)
#'
#' # DBSCAN with auto-detected eps
#' db <- cluster_tcrs(sub, "mouse", method = "dbscan")
#'
#' # Leiden on KNN graph
#' lei <- cluster_tcrs(sub, "mouse", method = "leiden", resolution = 1.0)
#'
#' # Precomputed distance matrix
#' dm <- tcrdist_matrix(sub, "mouse")
#' clusters <- cluster_tcrs(dist_matrix = dm, method = "kmedoids", k = 3)
#' }
#'
#' @seealso \code{\link{tcrdist_hclust}}, \code{\link{neighborhood_test}},
#'   \code{\link{compute_tcrdist_umap}}
#' @export
cluster_tcrs <- function(tcr_df = NULL, organism = NULL,
                          method = c("hierarchical", "leiden", "louvain",
                                     "dbscan", "kmedoids"),
                          dist_matrix = NULL,
                          k = NULL, h = NULL,
                          hclust_method = "average",
                          resolution = 1.0,
                          n_neighbors = 10L,
                          eps = NULL,
                          min_pts = 5L) {
    method <- match.arg(method)

    switch(method,
        hierarchical = {
            if (is.null(k) && is.null(h)) {
                stop("cluster_tcrs: exactly one of 'k' or 'h' must be ",
                     "specified for method='hierarchical'", call. = FALSE)
            }
            if (!is.null(k) && !is.null(h)) {
                stop("cluster_tcrs: specify either 'k' or 'h', not both",
                     call. = FALSE)
            }
            dm <- .get_dist_matrix(tcr_df, organism, dist_matrix)
            hc <- stats::hclust(stats::as.dist(dm), method = hclust_method)
            as.integer(stats::cutree(hc, k = k, h = h))
        },

        leiden = , louvain = {
            if (!requireNamespace("igraph", quietly = TRUE)) {
                stop("Package 'igraph' is required for method='", method,
                     "'. Install it with: install.packages(\"igraph\")",
                     call. = FALSE)
            }
            if (is.null(tcr_df) || is.null(organism)) {
                stop("cluster_tcrs: 'tcr_df' and 'organism' are required ",
                     "for method='", method, "'", call. = FALSE)
            }
            n <- nrow(tcr_df)
            knn_k <- min(as.integer(n_neighbors), n - 1L)
            knn <- tcrdist_knn(tcr_df, organism, K = knn_k)
            graph <- .build_knn_graph(knn$knn_indices, knn$knn_distances, n)
            clusters <- .run_clustering(graph, resolution = resolution,
                                         method = method)
            # .run_clustering returns 0-based; convert to 1-based
            as.integer(clusters + 1L)
        },

        dbscan = {
            if (!requireNamespace("dbscan", quietly = TRUE)) {
                stop("Package 'dbscan' is required for method='dbscan'. ",
                     "Install it with: install.packages(\"dbscan\")",
                     call. = FALSE)
            }
            dm <- .get_dist_matrix(tcr_df, organism, dist_matrix)
            if (is.null(eps)) {
                eps <- .dbscan_auto_eps(dm, min_pts)
            }
            res <- dbscan::dbscan(stats::as.dist(dm), eps = eps,
                                   minPts = min_pts)
            as.integer(res$cluster)
        },

        kmedoids = {
            if (!requireNamespace("cluster", quietly = TRUE)) {
                stop("Package 'cluster' is required for method='kmedoids'. ",
                     "Install it with: install.packages(\"cluster\")",
                     call. = FALSE)
            }
            if (is.null(k)) {
                stop("cluster_tcrs: 'k' is required for method='kmedoids'",
                     call. = FALSE)
            }
            dm <- .get_dist_matrix(tcr_df, organism, dist_matrix)
            res <- cluster::pam(stats::as.dist(dm), k = k, diss = TRUE)
            as.integer(res$clustering)
        }
    )
}


# ---------------------------------------------------------------------------
# neighborhood_test  (exported)
# ---------------------------------------------------------------------------

#' Test association between a variable and TCR neighborhoods
#'
#' For each TCR, tests whether a categorical variable is non-randomly
#' distributed among its TCRdist neighbors compared to the full repertoire.
#' Supports Fisher's exact test (binary variables) and chi-squared test
#' (multi-category variables).
#'
#' @details
#' For each TCR \eqn{i}, the test constructs a 2x2 (Fisher) or 2xK
#' (chi-squared) contingency table comparing category frequencies inside
#' the neighborhood (TCRs within \code{radius}) versus outside. The null
#' hypothesis is that the neighborhood is a random sample of the full
#' repertoire with respect to the variable.
#'
#' P-values are adjusted across all N tests using the method specified by
#' \code{p_adjust_method} (default: Benjamini-Hochberg, which controls
#' the false discovery rate). Note that TCR neighborhoods are spatially
#' correlated (nearby TCRs share neighbors), so the effective number of
#' independent tests is smaller than N. BH remains a reasonable choice
#' but may be conservative.
#'
#' @references
#' Dash, P. et al. (2017). Quantifiable predictive features define
#' epitope-specific T cell receptor repertoires. \emph{Nature}, 547, 89--93.
#'
#' @param tcr_df Data.frame with TCR columns (optional if \code{dist_matrix}
#'   is provided).
#' @param organism Character string (\code{"human"} or \code{"mouse"})
#'   (optional if \code{dist_matrix} is provided).
#' @param variable Character or factor vector of length \code{nrow(tcr_df)}.
#'   The categorical variable to test.
#' @param radius Numeric. Maximum TCRdist for neighborhood membership.
#'   Default \code{50}.
#' @param test Character string. \code{"fisher"} (default, for binary) or
#'   \code{"chisq"} (for multi-category).
#' @param p_adjust_method Character string. Method for \code{stats::p.adjust()}.
#'   Default \code{"BH"} (Benjamini-Hochberg).
#' @param dist_matrix Optional precomputed distance matrix. If provided,
#'   \code{tcr_df} and \code{organism} are not used for distance computation.
#'
#' @return A data.frame with one row per TCR and columns:
#'   \describe{
#'     \item{\code{index}}{Row index in \code{tcr_df}.}
#'     \item{\code{n_neighbors}}{Number of neighbors within radius.}
#'     \item{\code{p_value}}{Raw test p-value.}
#'     \item{\code{p_adjusted}}{Adjusted p-value.}
#'     \item{\code{odds_ratio}}{Odds ratio (Fisher only, NA for chi-sq).}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- neighborhood_test(tcr_df, "human",
#'                              variable = tcr_df$epitope, radius = 50)
#' significant <- result[result$p_adjusted < 0.05, ]
#' }
#'
#' @seealso \code{\link{tcrdist_hclust}}, \code{\link{cluster_tcrs}}, \code{\link{find_clumping}}
#' @export
neighborhood_test <- function(tcr_df = NULL, organism = NULL, variable,
                                radius = 50,
                                test = c("fisher", "chisq"),
                                p_adjust_method = "BH",
                                dist_matrix = NULL) {
    test <- match.arg(test)
    n <- length(variable)
    stopifnot(n >= 1L)

    variable <- as.factor(variable)
    lvls <- levels(variable)
    n_lvls <- length(lvls)

    if (test == "fisher" && n_lvls != 2L) {
        stop("neighborhood_test: Fisher's exact test requires exactly 2 ",
             "categories. Use test='chisq' for ", n_lvls, " categories.",
             call. = FALSE)
    }

    # Compute dense distance matrix
    dist_mat <- .get_dist_matrix(tcr_df, organism, dist_matrix)

    # Tally neighborhoods
    tally <- .neighborhood_tally(dist_mat, variable, radius)

    # Global category totals
    global_counts <- table(variable)

    # Run test for each TCR
    p_values <- numeric(n)
    odds_ratios <- rep(NA_real_, n)
    n_neighbors <- integer(n)

    for (i in seq_len(n)) {
        nbr_counts <- tally[i, ]
        n_nbr <- sum(nbr_counts)
        n_neighbors[i] <- n_nbr
        non_nbr_counts <- as.integer(global_counts) - nbr_counts

        if (n_nbr <= 1L || all(nbr_counts == 0L)) {
            p_values[i] <- 1.0
            next
        }

        if (test == "fisher") {
            # 2x2 table: category x neighborhood membership
            ct <- matrix(c(nbr_counts[1L], nbr_counts[2L],
                           non_nbr_counts[1L], non_nbr_counts[2L]),
                         nrow = 2, byrow = TRUE)
            ft <- stats::fisher.test(ct)
            p_values[i] <- ft$p.value
            odds_ratios[i] <- ft$estimate
        } else {
            ct <- rbind(nbr_counts, non_nbr_counts)
            # Suppress warnings for small expected counts
            ct_result <- suppressWarnings(stats::chisq.test(ct))
            p_values[i] <- ct_result$p.value
        }
    }

    p_adjusted <- stats::p.adjust(p_values, method = p_adjust_method)

    data.frame(
        index = seq_len(n),
        n_neighbors = n_neighbors,
        p_value = p_values,
        p_adjusted = p_adjusted,
        odds_ratio = odds_ratios,
        stringsAsFactors = FALSE
    )
}
