# TCR clone network construction and visualization
#
# Builds igraph networks from TCRdist distances and plots them with ggplot2.


# ---------------------------------------------------------------------------
# .find_distance_valley
# ---------------------------------------------------------------------------

#' Find the valley between the two major peaks of a bimodal distribution
#'
#' Uses kernel density estimation to locate the minimum between the two
#' dominant peaks of the typically bimodal TCRdist distribution.  First
#' identifies the two tallest local maxima, then finds the deepest local
#' minimum between them.  Falls back to a 2-component Gaussian mixture
#' model (EM) when KDE peak detection fails.
#'
#' @param distances Numeric vector of pairwise distances.
#' @param n_points Integer.  Number of density estimation points.
#' @param adjust Numeric.  Bandwidth adjustment for \code{density()}.
#' @return Numeric scalar.  The distance at the valley.
#' @keywords internal
.find_distance_valley <- function(distances, n_points = 512L, adjust = 1.0) {
    if (length(distances) < 10L) {
        warning(".find_distance_valley: too few distances, using median",
                call. = FALSE)
        return(stats::median(distances))
    }

    # --- KDE approach: find two tallest peaks, valley between them -----------
    d <- stats::density(distances, n = n_points, adjust = adjust)
    n <- length(d$y)

    # Local maxima
    is_max <- logical(n)
    for (i in 2L:(n - 1L)) {
        is_max[i] <- (d$y[i] > d$y[i - 1L]) && (d$y[i] > d$y[i + 1L])
    }
    max_idx <- which(is_max)

    # Local minima
    is_min <- logical(n)
    for (i in 2L:(n - 1L)) {
        is_min[i] <- (d$y[i] < d$y[i - 1L]) && (d$y[i] < d$y[i + 1L])
    }
    min_idx <- which(is_min)

    if (length(max_idx) >= 2L && length(min_idx) >= 1L) {
        # Sort peaks by height (tallest first)
        peak_order <- order(d$y[max_idx], decreasing = TRUE)
        top2 <- sort(max_idx[peak_order[1:2]])  # sorted by position

        # Find minima that lie between the two peaks
        between <- min_idx[min_idx > top2[1L] & min_idx < top2[2L]]

        if (length(between) > 0L) {
            # Deepest minimum between the two dominant peaks
            best <- between[which.min(d$y[between])]
            return(d$x[best])
        }
    }

    # --- GMM fallback: 2-component EM ----------------------------------------
    gmm <- .gmm_2_component(distances)
    if (!is.null(gmm)) return(gmm)

    # --- Final fallback: unimodal --------------------------------------------
    fallback <- stats::quantile(distances, 0.10, names = FALSE)
    message("No valley found in distance distribution (unimodal); ",
            "using 10th percentile: ", round(fallback, 1))
    fallback
}


#' Fit a 2-component Gaussian mixture via EM and find the crossover point
#'
#' Runs expectation-maximization for a mixture of two univariate Gaussians.
#' Returns the point where the posterior probabilities are equal (the
#' "crossover"), which corresponds to the valley between components.
#'
#' @param x Numeric vector of distances.
#' @param max_iter Integer.  Maximum EM iterations.
#' @param tol Numeric.  Convergence tolerance on log-likelihood.
#' @return Numeric scalar (crossover point) or \code{NULL} on failure.
#' @keywords internal
.gmm_2_component <- function(x, max_iter = 200L, tol = 1e-6) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 20L) return(NULL)

    # Initialise with K-means-style split at the median
    med <- stats::median(x)
    idx1 <- x <= med
    mu1  <- mean(x[idx1]);  mu2  <- mean(x[!idx1])
    sd1  <- stats::sd(x[idx1]);  sd2  <- stats::sd(x[!idx1])
    if (sd1 < 1e-10) sd1 <- 1;  if (sd2 < 1e-10) sd2 <- 1
    pi1  <- sum(idx1) / n;  pi2  <- 1 - pi1

    prev_ll <- -Inf
    for (iter in seq_len(max_iter)) {
        # E-step
        d1 <- pi1 * stats::dnorm(x, mu1, sd1)
        d2 <- pi2 * stats::dnorm(x, mu2, sd2)
        total <- d1 + d2
        total[total < 1e-300] <- 1e-300
        gamma1 <- d1 / total

        # M-step
        n1  <- sum(gamma1);  n2  <- n - n1
        if (n1 < 2 || n2 < 2) return(NULL)
        pi1 <- n1 / n;  pi2 <- 1 - pi1
        mu1 <- sum(gamma1 * x) / n1
        mu2 <- sum((1 - gamma1) * x) / n2
        sd1 <- sqrt(sum(gamma1 * (x - mu1)^2) / n1)
        sd2 <- sqrt(sum((1 - gamma1) * (x - mu2)^2) / n2)
        if (sd1 < 1e-10) sd1 <- 1e-10
        if (sd2 < 1e-10) sd2 <- 1e-10

        ll <- sum(log(pi1 * stats::dnorm(x, mu1, sd1) +
                      pi2 * stats::dnorm(x, mu2, sd2)))
        if (is.finite(ll) && abs(ll - prev_ll) < tol) break
        prev_ll <- ll
    }

    # Ensure mu1 < mu2
    if (mu1 > mu2) {
        tmp <- mu1; mu1 <- mu2; mu2 <- tmp
        tmp <- sd1; sd1 <- sd2; sd2 <- tmp
        tmp <- pi1; pi1 <- pi2; pi2 <- tmp
    }

    # Find crossover: where pi1*N(x|mu1,sd1) == pi2*N(x|mu2,sd2)
    # Search on a fine grid between the two means
    grid <- seq(mu1, mu2, length.out = 1000L)
    diff <- abs(pi1 * stats::dnorm(grid, mu1, sd1) -
                pi2 * stats::dnorm(grid, mu2, sd2))
    crossover <- grid[which.min(diff)]
    crossover
}


# ---------------------------------------------------------------------------
# .sparse_to_edgelist
# ---------------------------------------------------------------------------

#' Extract edge list from a sparse distance matrix
#'
#' Converts upper triangle of a symmetric sparse matrix to an edge list
#' with distance values.
#'
#' @param sp A \code{dgCMatrix} (symmetric sparse distance matrix).
#' @return A \code{data.frame} with columns \code{from}, \code{to},
#'   \code{distance} (all entries from upper triangle with dist > 0).
#' @keywords internal
.sparse_to_edgelist <- function(sp) {
    sp_t <- methods::as(sp, "TsparseMatrix")
    # 0-based indices; take upper triangle only
    upper <- sp_t@i < sp_t@j
    data.frame(
        from     = sp_t@i[upper] + 1L,
        to       = sp_t@j[upper] + 1L,
        distance = sp_t@x[upper],
        stringsAsFactors = FALSE
    )
}


# ---------------------------------------------------------------------------
# .jitter_overlapping
# ---------------------------------------------------------------------------

#' Jitter overlapping layout coordinates
#'
#' Finds groups of nodes that share identical (x, y) positions and
#' spreads each group into a small circle.  The radius is proportional
#' to the overall layout extent so the offset is visible but not
#' disruptive.
#'
#' @param coords Numeric matrix (N x 2) of layout coordinates.
#' @param frac Numeric.  Jitter radius as a fraction of the layout
#'   extent.  Default \code{0.015}.
#' @return Numeric matrix (N x 2) with adjusted coordinates.
#' @keywords internal
.jitter_overlapping <- function(coords, frac = 0.015) {
    key <- paste(coords[, 1L], coords[, 2L], sep = "\x01")
    groups <- split(seq_len(nrow(coords)), key)

    # Only process groups with 2+ overlapping nodes
    dups <- groups[lengths(groups) >= 2L]
    if (length(dups) == 0L) return(coords)

    # Scale radius to the layout extent
    x_range <- diff(range(coords[, 1L]))
    y_range <- diff(range(coords[, 2L]))
    extent  <- max(x_range, y_range, 1e-6)
    radius  <- frac * extent

    for (idx in dups) {
        k <- length(idx)
        angles <- seq(0, 2 * pi, length.out = k + 1L)[seq_len(k)]
        coords[idx, 1L] <- coords[idx, 1L] + radius * cos(angles)
        coords[idx, 2L] <- coords[idx, 2L] + radius * sin(angles)
    }
    coords
}


# ===========================================================================
# compute_tcr_network
# ===========================================================================

#' Build a TCR clone network from distance data
#'
#' Constructs an \pkg{igraph} graph where vertices are TCR clones and edges
#' connect clones within a distance threshold.  Edge weights encode similarity
#' via exponential decay: \code{exp(-distance / scale)}.
#'
#' When \code{threshold = NULL}, the threshold is auto-detected by finding the
#' valley between the two peaks of the (typically bimodal) TCRdist
#' distribution, using kernel density estimation on a subsample.
#'
#' A distance distribution plot with the threshold marked is automatically
#' displayed when \pkg{ggplot2} is available.  The plot is also stored in the
#' returned list as \code{dist_plot}.
#'
#' @param tcrs Data.frame with TCR columns (\code{va}, \code{cdr3a},
#'   \code{vb}, \code{cdr3b}).  All columns are attached as vertex attributes
#'   on the output graph.
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#'   Ignored when \code{dist_matrix} is provided.
#' @param threshold Numeric or \code{NULL}.  Maximum TCRdist for an edge.
#'   If \code{NULL}, auto-detected from the distance distribution valley.
#' @param dist_matrix Optional precomputed N x N distance matrix (e.g. from
#'   \code{tcrdist_matrix()} or \code{TCRrep@@paired_dist}).  If provided,
#'   \code{organism} is ignored and distances are taken from this matrix.
#' @param scale Numeric or \code{NULL}.  Scale parameter for the similarity
#'   transform \code{exp(-dist / scale)}.  If \code{NULL}, defaults to
#'   \code{threshold / 4}.
#' @param min_edges Integer.  Minimum number of edges a vertex must have to
#'   be kept.  \code{0} (default) keeps all vertices.  \code{1} removes
#'   singletons (isolated nodes), \code{2} removes vertices with fewer than
#'   2 edges, etc.
#' @param jitter Logical.  If \code{TRUE} (default), slightly offset
#'   overlapping nodes (e.g. identical clones with distance 0) so they are
#'   individually visible instead of stacking on top of each other.
#' @param layout Character.  Layout algorithm: \code{"fr"}
#'   (Fruchterman-Reingold), \code{"kk"} (Kamada-Kawai), \code{"drl"},
#'   \code{"circle"}, or \code{"grid"}.  Default \code{"fr"}.
#' @param seed Integer or \code{NULL}.  Random seed for layout
#'   reproducibility.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{graph}}{An \code{igraph} object.  Vertex attributes include
#'       all columns from \code{tcrs}.  Edge attributes: \code{weight}
#'       (similarity) and \code{distance} (TCRdist).}
#'     \item{\code{layout}}{Numeric matrix (N x 2) of layout coordinates.}
#'     \item{\code{threshold}}{Numeric.  Threshold used.}
#'     \item{\code{scale}}{Numeric.  Scale parameter used.}
#'     \item{\code{n_components}}{Integer.  Number of connected components.}
#'     \item{\code{n_edges}}{Integer.  Number of edges in the graph.}
#'     \item{\code{dist_plot}}{\code{ggplot} object showing the distance
#'       distribution with the threshold line, or \code{NULL} if
#'       \pkg{ggplot2} is not available.}
#'   }
#'
#' @examples
#' \donttest{
#' data(dash)
#' sub <- dash[1:200, ]
#' net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42)
#' net$n_components
#' net$n_edges
#' }
#'
#' @seealso \code{\link{plot_tcr_network}}, \code{\link{tcrdist_sparse}},
#'   \code{\link{tcrdist_matrix}}
#' @export
compute_tcr_network <- function(tcrs,
                                 organism = NULL,
                                 threshold = NULL,
                                 dist_matrix = NULL,
                                 scale = NULL,
                                 min_edges = 0L,
                                 jitter = TRUE,
                                 layout = "fr",
                                 seed = NULL) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required for compute_tcr_network(). ",
             "Install it with: install.packages(\"igraph\")",
             call. = FALSE)
    }
    if (!is.data.frame(tcrs)) {
        stop("compute_tcr_network: 'tcrs' must be a data.frame", call. = FALSE)
    }
    n <- nrow(tcrs)
    if (n < 2L) {
        stop("compute_tcr_network: need at least 2 clones", call. = FALSE)
    }

    # ---- Get edges -----------------------------------------------------------
    dists_vec <- NULL  # for distribution plot

    if (!is.null(dist_matrix)) {
        # From precomputed distance matrix
        dist_matrix <- as.matrix(dist_matrix)
        if (nrow(dist_matrix) != n || ncol(dist_matrix) != n) {
            stop("compute_tcr_network: 'dist_matrix' dimensions must match ",
                 "nrow(tcrs)", call. = FALSE)
        }

        dists_vec <- dist_matrix[upper.tri(dist_matrix)]

        if (is.null(threshold)) {
            threshold <- .find_distance_valley(dists_vec)
            message("Auto-detected distance threshold: ", round(threshold, 1))
        }

        # Extract edges from upper triangle within threshold
        idx <- which(dist_matrix <= threshold & upper.tri(dist_matrix),
                     arr.ind = TRUE)
        if (nrow(idx) > 0L) {
            edge_from <- idx[, 1L]
            edge_to   <- idx[, 2L]
            edge_dist <- dist_matrix[idx]
        } else {
            edge_from <- integer(0L)
            edge_to   <- integer(0L)
            edge_dist <- numeric(0L)
        }
    } else {
        # Compute sparse distances
        if (is.null(organism)) {
            stop("compute_tcr_network: 'organism' is required when ",
                 "'dist_matrix' is not provided", call. = FALSE)
        }

        # Subsample to get distance distribution (for threshold + plot)
        n_sub <- min(500L, n)
        sub_idx <- if (n_sub < n) sample(n, n_sub) else seq_len(n)
        sub_dists <- tcrdist_matrix(tcrs[sub_idx, ], organism)
        dists_vec <- sub_dists[upper.tri(sub_dists)]

        if (is.null(threshold)) {
            threshold <- .find_distance_valley(dists_vec)
            message("Auto-detected distance threshold: ", round(threshold, 1))
        }

        sp <- tcrdist_sparse(tcrs, organism, threshold = threshold)
        el <- .sparse_to_edgelist(sp)
        edge_from <- el$from
        edge_to   <- el$to
        edge_dist <- el$distance

        # Handle zero-distance pairs (sparseMatrix drops explicit zeros)
        chain_keys <- paste(tcrs$va, tcrs$cdr3a, tcrs$vb, tcrs$cdr3b,
                            sep = "\x01")
        dup_groups <- which(duplicated(chain_keys) | duplicated(chain_keys,
                                                                 fromLast = TRUE))
        if (length(dup_groups) > 0L) {
            key_to_idx <- split(seq_len(n), chain_keys)
            for (group in key_to_idx) {
                if (length(group) < 2L) next
                pairs <- utils::combn(sort(group), 2L)
                edge_from <- c(edge_from, pairs[1L, ])
                edge_to   <- c(edge_to,   pairs[2L, ])
                edge_dist <- c(edge_dist, rep(0, ncol(pairs)))
            }
        }
    }

    # ---- Build igraph --------------------------------------------------------
    if (is.null(scale)) scale <- threshold / 4
    if (scale <= 0) scale <- 1  # safety

    g <- igraph::make_empty_graph(n = n, directed = FALSE)

    if (length(edge_from) > 0L) {
        edge_weight <- exp(-edge_dist / scale)
        edge_pairs  <- rbind(edge_from, edge_to)
        g <- igraph::add_edges(g, as.vector(edge_pairs))
        igraph::E(g)$weight   <- edge_weight
        igraph::E(g)$distance <- edge_dist
    }

    # Attach clone metadata as vertex attributes
    for (col in colnames(tcrs)) {
        g <- igraph::set_vertex_attr(g, col, value = tcrs[[col]])
    }

    # ---- Prune low-degree vertices ------------------------------------------
    min_edges <- as.integer(min_edges)
    if (min_edges > 0L) {
        repeat {
            deg <- igraph::degree(g)
            drop <- which(deg < min_edges)
            if (length(drop) == 0L) break
            keep <- which(deg >= min_edges)
            if (length(keep) == 0L) {
                warning("compute_tcr_network: min_edges=", min_edges,
                        " removed all vertices; returning unpruned graph",
                        call. = FALSE)
                break
            }
            g <- igraph::induced_subgraph(g, keep)
        }
    }

    # ---- Layout --------------------------------------------------------------
    if (!is.null(seed)) set.seed(seed)

    layout_fn <- switch(
        layout,
        "fr"     = igraph::layout_with_fr,
        "kk"     = igraph::layout_with_kk,
        "drl"    = igraph::layout_with_drl,
        "circle" = igraph::layout_in_circle,
        "grid"   = igraph::layout_on_grid,
        stop(sprintf(
            "Unknown layout: '%s'. Use one of: fr, kk, drl, circle, grid",
            layout
        ), call. = FALSE)
    )

    # FR/KK use weights for spring strength; others ignore them
    coords <- if (layout %in% c("fr", "kk") && igraph::ecount(g) > 0L) {
        layout_fn(g, weights = igraph::E(g)$weight)
    } else {
        layout_fn(g)
    }

    # ---- Jitter overlapping nodes -------------------------------------------
    if (jitter) {
        coords <- .jitter_overlapping(coords)
    }

    # ---- Distance distribution plot with threshold line ---------------------
    dist_plot <- NULL
    if (!is.null(dists_vec) && length(dists_vec) > 0L &&
        requireNamespace("ggplot2", quietly = TRUE)) {
        dist_plot <- plot_distance_distribution(
            dists_vec, threshold = threshold,
            title = "TCRdist distribution with network threshold"
        )
        print(dist_plot)
    }

    list(
        graph        = g,
        layout       = coords,
        threshold    = threshold,
        scale        = scale,
        n_components = igraph::components(g)$no,
        n_edges      = igraph::ecount(g),
        dist_plot    = dist_plot
    )
}


# ===========================================================================
# plot_tcr_network
# ===========================================================================

#' Plot a TCR clone network
#'
#' Draws a network graph of TCR clones with ggplot2, using layout coordinates
#' from \code{\link{compute_tcr_network}}.  Vertices are colored by a metadata
#' variable (e.g. epitope, cluster) and edges are drawn as semi-transparent
#' segments.
#'
#' @param network List returned by \code{\link{compute_tcr_network}}.
#' @param color_by Character string naming a vertex attribute, OR a vector
#'   of length N.  If a single string matching a vertex attribute name, that
#'   attribute is used.  If \code{NULL}, all vertices are gray.
#' @param vertex_size Numeric.  Point size.  Default \code{3}.
#' @param edge_alpha Numeric.  Edge transparency.  Default \code{0.15}.
#' @param edge_width Numeric.  Edge line width.  Default \code{0.3}.
#' @param title Optional plot title.
#' @param palette Character vector of colors to override the default palette.
#' @param legend_title Optional legend title.
#' @param show_labels Logical.  If \code{TRUE} and coloring is categorical,
#'   add centroid labels.  Default \code{FALSE}.
#' @param label_size Numeric.  Label text size.  Default \code{3}.
#' @param na_color Character.  Color for NA values.  Default
#'   \code{"#DDDDDD"}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' data(dash)
#' sub <- dash[1:200, ]
#' net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42)
#' plot_tcr_network(net, color_by = "epitope", title = "TCR network")
#' }
#'
#' @seealso \code{\link{compute_tcr_network}}, \code{\link{plot_tcr_scatter}}
#' @export
plot_tcr_network <- function(network,
                              color_by = NULL,
                              vertex_size = 3,
                              edge_alpha = 0.15,
                              edge_width = 0.3,
                              title = NULL,
                              palette = NULL,
                              legend_title = NULL,
                              show_labels = FALSE,
                              label_size = 3,
                              na_color = "#DDDDDD") {
    .check_ggplot2("plot_tcr_network()")

    g      <- network$graph
    coords <- network$layout
    n      <- igraph::vcount(g)

    # ---- Vertex data ---------------------------------------------------------
    vertex_df <- data.frame(
        x = coords[, 1L],
        y = coords[, 2L],
        stringsAsFactors = FALSE
    )

    # ---- Resolve color_by ----------------------------------------------------
    if (!is.null(color_by)) {
        if (is.character(color_by) && length(color_by) == 1L &&
            color_by %in% igraph::vertex_attr_names(g)) {
            color_vals <- igraph::vertex_attr(g, color_by)
            if (is.null(legend_title)) legend_title <- color_by
        } else {
            color_vals <- color_by
        }
        stopifnot(length(color_vals) == n)
        vertex_df$color <- color_vals
    }

    # ---- Edge data -----------------------------------------------------------
    edge_df <- NULL
    if (igraph::ecount(g) > 0L) {
        el <- igraph::as_edgelist(g, names = FALSE)
        edge_df <- data.frame(
            x    = coords[el[, 1L], 1L],
            y    = coords[el[, 1L], 2L],
            xend = coords[el[, 2L], 1L],
            yend = coords[el[, 2L], 2L],
            stringsAsFactors = FALSE
        )
    }

    # ---- Build plot ----------------------------------------------------------
    p <- ggplot2::ggplot()

    # Edges
    if (!is.null(edge_df) && nrow(edge_df) > 0L) {
        p <- p + ggplot2::geom_segment(
            data = edge_df,
            ggplot2::aes(x = .data$x, y = .data$y,
                         xend = .data$xend, yend = .data$yend),
            alpha = edge_alpha, linewidth = edge_width, color = "gray70"
        )
    }

    # Vertices
    if (is.null(color_by)) {
        p <- p + ggplot2::geom_point(
            data = vertex_df,
            ggplot2::aes(x = .data$x, y = .data$y),
            size = vertex_size, color = "gray50"
        )
    } else {
        is_categorical <- is.factor(vertex_df$color) ||
            is.character(vertex_df$color)

        p <- p + ggplot2::geom_point(
            data = vertex_df,
            ggplot2::aes(x = .data$x, y = .data$y, color = .data$color),
            size = vertex_size
        )

        if (is_categorical) {
            lvls <- if (is.factor(vertex_df$color)) levels(vertex_df$color)
                    else sort(unique(vertex_df$color))
            if (is.null(palette)) {
                pal <- .tcrdistR_palette(length(lvls))
            } else {
                pal <- rep_len(palette, length(lvls))
            }
            names(pal) <- lvls
            p <- p + ggplot2::scale_color_manual(values = pal,
                                                  na.value = na_color)

            if (show_labels) {
                centroids <- stats::aggregate(
                    cbind(x, y) ~ color, data = vertex_df, FUN = mean)
                if (requireNamespace("ggrepel", quietly = TRUE)) {
                    p <- p + ggrepel::geom_label_repel(
                        data = centroids,
                        ggplot2::aes(x = .data$x, y = .data$y,
                                     label = .data$color),
                        size = label_size, fontface = "bold",
                        fill = "white", alpha = 0.8,
                        label.size = 0.2, label.padding = 0.15,
                        max.overlaps = Inf, show.legend = FALSE)
                } else {
                    p <- p + ggplot2::geom_label(
                        data = centroids,
                        ggplot2::aes(x = .data$x, y = .data$y,
                                     label = .data$color),
                        size = label_size, fontface = "bold",
                        fill = "white", alpha = 0.8,
                        label.size = 0.2,
                        label.padding = ggplot2::unit(0.15, "lines"),
                        show.legend = FALSE)
                }
            }
        } else {
            p <- p + ggplot2::scale_color_viridis_c(na.value = na_color)
        }
    }

    p + ggplot2::labs(title = title, color = legend_title) +
        ggplot2::coord_fixed() +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text  = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        )
}
