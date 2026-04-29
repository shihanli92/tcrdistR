# Generic 2D scatter plot for PCA/UMAP embeddings
#
# Ported from rconga R/plot_core.R plot_conga_umap().
# Dependencies: ggplot2 (Suggests), ggrepel (Suggests, optional)


# ---------------------------------------------------------------------------
# plot_tcr_scatter  (exported)
# ---------------------------------------------------------------------------

#' Plot a 2D scatter of TCR embeddings
#'
#' Generic scatter plot for kernel PCA, UMAP, or any 2D embedding of TCR
#' repertoire data. Supports continuous (viridis) and categorical
#' (tab10/tab20) coloring with optional centroid labels.
#'
#' @param coords Numeric matrix with 2 columns (embedding coordinates).
#' @param color_by Optional vector of length \code{nrow(coords)}. Numeric
#'   for continuous coloring (viridis), factor/character for categorical
#'   coloring (tab10/tab20). \code{NULL} plots all points in gray.
#' @param title Optional plot title.
#' @param point_size Numeric. Point size. Default \code{1}.
#' @param alpha Numeric. Point opacity. Default \code{1}.
#' @param axis_label_prefix Character string. Prefix for axis labels.
#'   Default \code{"KPCA"}.
#' @param legend_title Optional legend title.
#' @param palette Character vector of colors to override default palette.
#' @param show_labels Logical. If \code{TRUE} and \code{color_by} is
#'   categorical, add centroid labels. Default \code{FALSE}.
#' @param label_size Numeric. Label text size. Default \code{3}.
#' @param na_color Character. Color for NA values. Default \code{"#DDDDDD"}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' kpca <- compute_tcrdist_kernel_pca(tcr_df, "human")
#' plot_tcr_scatter(kpca$embeddings[, 1:2], color_by = tcr_df$epitope)
#' }
#'
#' @export
plot_tcr_scatter <- function(coords, color_by = NULL, title = NULL,
                              point_size = 1, alpha = 1,
                              axis_label_prefix = "KPCA",
                              legend_title = NULL,
                              palette = NULL,
                              show_labels = FALSE,
                              label_size = 3,
                              na_color = "#DDDDDD") {
    .check_ggplot2("plot_tcr_scatter()")

    coords <- as.matrix(coords)
    stopifnot(ncol(coords) == 2L)

    df <- data.frame(
        x = coords[, 1L],
        y = coords[, 2L],
        stringsAsFactors = FALSE
    )

    if (is.null(color_by)) {
        # All gray
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
            ggplot2::geom_point(size = point_size, alpha = alpha,
                                color = "gray50")
    } else {
        stopifnot(length(color_by) == nrow(coords))
        is_categorical <- is.factor(color_by) || is.character(color_by)
        df$color <- color_by

        # Sort so high-value points draw on top (continuous only)
        if (!is_categorical) {
            df <- df[order(df$color, na.last = FALSE), ]
        }

        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                                color = .data$color)) +
            ggplot2::geom_point(size = point_size, alpha = alpha)

        if (is_categorical) {
            lvls <- if (is.factor(color_by)) levels(color_by)
                    else sort(unique(color_by))
            if (is.null(palette)) {
                pal <- .tcrdistR_palette(length(lvls))
            } else {
                pal <- rep_len(palette, length(lvls))
            }
            names(pal) <- lvls
            p <- p + ggplot2::scale_color_manual(values = pal,
                                                  na.value = na_color)

            # Centroid labels
            if (show_labels) {
                centroids <- stats::aggregate(
                    cbind(x, y) ~ color, data = df, FUN = mean)
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

    p + ggplot2::labs(
            x = paste0(axis_label_prefix, "1"),
            y = paste0(axis_label_prefix, "2"),
            title = title,
            color = legend_title
        ) +
        ggplot2::coord_fixed() +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        )
}
