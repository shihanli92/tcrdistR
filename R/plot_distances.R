# Distance matrix visualization: heatmaps, dendrograms, distributions
#
# Ported from rconga R/plot_core.R (Python CoNGA plotting.py).
# Dependencies: ggplot2 (Suggests)


# ---------------------------------------------------------------------------
# .hclust_to_segments  (internal)
# ---------------------------------------------------------------------------

#' Convert an hclust object to ggplot2 segment coordinates
#'
#' Recursively traverses a dendrogram and extracts (x, y, xend, yend)
#' segment data suitable for \code{ggplot2::geom_segment()}.
#'
#' @param hc An \code{hclust} object.
#' @return A data.frame with columns \code{x}, \code{y}, \code{xend},
#'   \code{yend}.
#' @keywords internal
#' @noRd
.hclust_to_segments <- function(hc) {
    dd <- stats::as.dendrogram(hc)

    # Label leaves with their plot positions
    dd <- stats::dendrapply(dd, function(node) {
        if (stats::is.leaf(node)) {
            leaf_val <- as.integer(node)
            pos <- match(leaf_val, hc$order)
            attr(node, "label_x") <- pos
        }
        node
    })

    seg_env <- new.env(parent = emptyenv())
    seg_env$segs <- list()

    .traverse <- function(node, segments_env) {
        if (stats::is.leaf(node)) {
            return(list(
                x    = as.numeric(attr(node, "label_x") %||%
                    match(as.integer(node), hc$order)),
                ybot = 0
            ))
        }

        children <- lapply(seq_along(node), function(i) {
            .traverse(node[[i]], segments_env)
        })

        h <- attr(node, "height")
        xs <- vapply(children, `[[`, numeric(1), "x")
        ybots <- vapply(children, `[[`, numeric(1), "ybot")

        # Vertical segments from each child up to merge height
        for (i in seq_along(children)) {
            segments_env$segs <- c(segments_env$segs, list(data.frame(
                x = xs[i], y = ybots[i], xend = xs[i], yend = h)))
        }
        # Horizontal segment connecting children
        segments_env$segs <- c(segments_env$segs, list(data.frame(
            x = min(xs), y = h, xend = max(xs), yend = h)))

        list(x = mean(xs), ybot = h)
    }

    .traverse(dd, seg_env)

    if (length(seg_env$segs) == 0L) {
        return(data.frame(x = numeric(0), y = numeric(0),
                          xend = numeric(0), yend = numeric(0)))
    }

    do.call(rbind, seg_env$segs)
}


# ---------------------------------------------------------------------------
# plot_tcrdist_heatmap  (exported)
# ---------------------------------------------------------------------------

#' Plot a TCRdist pairwise distance heatmap
#'
#' Renders a heatmap of a pairwise distance matrix using
#' \code{ggplot2::geom_tile()}. Optionally clusters rows and columns via
#' hierarchical clustering (UPGMA).
#'
#' @param dist_matrix Numeric matrix. Square pairwise distance matrix.
#' @param labels Character vector. Row/column labels. If \code{NULL},
#'   uses \code{rownames(dist_matrix)} or integer indices.
#' @param cluster Logical. If \code{TRUE} (default), reorder rows and columns
#'   by hierarchical clustering.
#' @param title Optional plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(c(0,10,20, 10,0,15, 20,15,0), nrow = 3)
#' plot_tcrdist_heatmap(mat, labels = c("TCR1", "TCR2", "TCR3"))
#' }
#'
#' @seealso \code{\link{plot_tcrdist_dendrogram}}, \code{\link{plot_distance_distribution}}, \code{\link{tcrdist_matrix}}
#' @export
plot_tcrdist_heatmap <- function(dist_matrix, labels = NULL,
                                  cluster = TRUE, title = NULL) {
    .check_ggplot2("plot_tcrdist_heatmap()")

    stopifnot(is.matrix(dist_matrix), nrow(dist_matrix) == ncol(dist_matrix))
    n <- nrow(dist_matrix)

    if (is.null(labels)) {
        labels <- rownames(dist_matrix) %||% as.character(seq_len(n))
    }

    # Cluster ordering
    if (cluster && n >= 2L) {
        hc <- stats::hclust(stats::as.dist(dist_matrix), method = "average")
        ord <- hc$order
    } else {
        ord <- seq_len(n)
    }

    labels_ordered <- labels[ord]

    # Build long-form data.frame
    df_list <- vector("list", n * n)
    idx <- 0L
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            idx <- idx + 1L
            ri <- ord[i]
            ci <- ord[j]
            df_list[[idx]] <- data.frame(
                row = i,
                col = j,
                row_label = labels_ordered[i],
                col_label = labels_ordered[j],
                distance = dist_matrix[ri, ci],
                stringsAsFactors = FALSE
            )
        }
    }
    df <- do.call(rbind, df_list)

    df$row_label <- factor(df$row_label, levels = rev(labels_ordered))
    df$col_label <- factor(df$col_label, levels = labels_ordered)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$col_label,
                                      y = .data$row_label,
                                      fill = .data$distance)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
            low = "#1f77b4", mid = "white", high = "#d62728",
            midpoint = stats::median(dist_matrix[upper.tri(dist_matrix)]),
            name = "Distance"
        ) +
        ggplot2::labs(x = NULL, y = NULL, title = title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                 size = 7),
            axis.text.y = ggplot2::element_text(size = 7),
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
        )
}


# ---------------------------------------------------------------------------
# plot_tcrdist_dendrogram  (exported)
# ---------------------------------------------------------------------------

#' Plot a TCRdist-based dendrogram
#'
#' Computes pairwise TCRdist distances, performs hierarchical clustering, and
#' renders a dendrogram with optional leaf coloring.
#'
#' @param tcr_df Data.frame with TCR columns (\code{va}, \code{vb},
#'   \code{cdr3a}, \code{cdr3b}).
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#' @param color_by Optional vector for leaf coloring. Numeric gives a viridis
#'   gradient; character/factor gives a qualitative palette. \code{NULL}
#'   colors all leaves gray.
#' @param method Clustering method for \code{stats::hclust()}. Default
#'   \code{"average"} (UPGMA).
#' @param max_tcrs Integer. If \code{nrow(tcr_df)} exceeds this, subsample.
#'   Default \code{500L}.
#' @param title Optional plot title.
#' @param point_size Numeric. Leaf point size. Default \code{2}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plot_tcrdist_dendrogram(tcr_df, "human", color_by = tcr_df$cluster)
#' }
#'
#' @seealso \code{\link{plot_tcrdist_heatmap}}, \code{\link{tcrdist_hclust}}, \code{\link{cluster_tcrs}}
#' @export
plot_tcrdist_dendrogram <- function(tcr_df, organism, color_by = NULL,
                                     method = "average", max_tcrs = 500L,
                                     title = NULL, point_size = 2) {
    .check_ggplot2("plot_tcrdist_dendrogram()")

    n <- nrow(tcr_df)
    if (n < 2L) {
        stop("Need at least 2 TCRs for a dendrogram.", call. = FALSE)
    }

    # Subsample if needed
    if (n > max_tcrs) {
        sel <- sample(n, max_tcrs)
        tcr_df <- tcr_df[sel, , drop = FALSE]
        if (!is.null(color_by)) color_by <- color_by[sel]
        n <- max_tcrs
    }

    dist_mat <- tcrdist_matrix(tcr_df, organism)
    hc <- stats::hclust(stats::as.dist(dist_mat), method = method)
    seg <- .hclust_to_segments(hc)

    leaf_order <- hc$order
    leaf_df <- data.frame(
        x = seq_along(leaf_order),
        y = 0,
        stringsAsFactors = FALSE
    )

    is_categorical <- FALSE
    if (!is.null(color_by)) {
        leaf_df$color <- color_by[leaf_order]
        is_categorical <- is.factor(color_by) || is.character(color_by)
        if (is_categorical) {
            leaf_df$color <- as.factor(leaf_df$color)
        }
    }

    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = seg,
            ggplot2::aes(x = .data$x, y = .data$y,
                         xend = .data$xend, yend = .data$yend),
            color = "gray40", linewidth = 0.3
        )

    if (!is.null(color_by)) {
        p <- p + ggplot2::geom_point(
            data = leaf_df,
            ggplot2::aes(x = .data$x, y = .data$y, color = .data$color),
            size = point_size
        )
        if (is_categorical) {
            lvls <- levels(leaf_df$color)
            pal <- .tcrdistR_palette(length(lvls))
            names(pal) <- lvls
            p <- p + ggplot2::scale_color_manual(values = pal)
        } else {
            p <- p + ggplot2::scale_color_viridis_c()
        }
    } else {
        p <- p + ggplot2::geom_point(
            data = leaf_df,
            ggplot2::aes(x = .data$x, y = .data$y),
            color = "gray50", size = point_size
        )
    }

    p + ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.05, 0.1))
        ) +
        ggplot2::labs(title = title, color = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        )
}


# ---------------------------------------------------------------------------
# plot_distance_distribution  (exported)
# ---------------------------------------------------------------------------

#' Plot the distribution of pairwise TCRdist distances
#'
#' Histogram with overlaid density curve of the upper triangle of a pairwise
#' distance matrix.  A vertical dashed line marks the median distance.
#' When \code{threshold} is provided, a red vertical line is drawn at that
#' value with a label.
#'
#' @param dist_matrix Numeric matrix or numeric vector.  If a square matrix,
#'   the upper triangle is extracted; if a vector, used directly.
#' @param title Optional plot title.
#' @param binwidth Numeric. Histogram bin width. If \code{NULL}, uses
#'   ggplot2 default.
#' @param threshold Numeric or \code{NULL}.  If not \code{NULL}, a vertical
#'   line is drawn at this distance value and labeled.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(c(0,10,20, 10,0,15, 20,15,0), nrow = 3)
#' plot_distance_distribution(mat, threshold = 12)
#' }
#'
#' @seealso \code{\link{plot_tcrdist_heatmap}}, \code{\link{tcrdist_matrix}},
#'   \code{\link{compute_tcr_network}}
#' @export
plot_distance_distribution <- function(dist_matrix, title = NULL,
                                        binwidth = NULL,
                                        threshold = NULL) {
    .check_ggplot2("plot_distance_distribution()")

    # Accept a vector or a square matrix
    if (is.matrix(dist_matrix)) {
        stopifnot(nrow(dist_matrix) == ncol(dist_matrix))
        dists <- dist_matrix[upper.tri(dist_matrix)]
    } else {
        dists <- as.numeric(dist_matrix)
    }

    if (length(dists) == 0L) {
        return(ggplot2::ggplot() +
                   ggplot2::theme_void() +
                   ggplot2::ggtitle(title %||% "No pairwise distances"))
    }

    df <- data.frame(distance = dists)
    med <- stats::median(dists)

    if (is.null(title)) {
        title <- "Pairwise distance distribution"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$distance))

    if (!is.null(binwidth)) {
        p <- p + ggplot2::geom_histogram(
            ggplot2::aes(y = ggplot2::after_stat(.data$density)),
            binwidth = binwidth, fill = "#1f77b4", alpha = 0.6)
    } else {
        p <- p + ggplot2::geom_histogram(
            ggplot2::aes(y = ggplot2::after_stat(.data$density)),
            bins = 30L, fill = "#1f77b4", alpha = 0.6)
    }

    p <- p + ggplot2::geom_density(color = "#d62728", linewidth = 0.8) +
        ggplot2::geom_vline(xintercept = med, linetype = "dashed",
                            color = "gray40")

    # Threshold line
    if (!is.null(threshold)) {
        p <- p + ggplot2::geom_vline(
            xintercept = threshold, color = "#d62728",
            linetype = "solid", linewidth = 0.9
        ) +
        ggplot2::annotate(
            "text", x = threshold, y = Inf,
            label = paste0("threshold = ", round(threshold, 1)),
            vjust = 2, hjust = -0.05, color = "#d62728",
            fontface = "bold", size = 3.5
        )
    }

    p + ggplot2::labs(x = "TCRdist", y = "Density", title = title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
        )
}
