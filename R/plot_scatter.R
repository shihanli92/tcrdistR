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
#' (tab10/tab20) coloring with optional centroid labels, faceting, and
#' clone highlighting.
#'
#' @param coords Numeric matrix with 2 columns (embedding coordinates).
#' @param color_by Optional vector of length \code{nrow(coords)}. Numeric
#'   for continuous coloring (viridis), factor/character for categorical
#'   coloring (tab10/tab20). \code{NULL} plots all points in gray. When
#'   \code{metadata} is provided, can be a column name string.
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
#' @param metadata Optional data.frame with \code{nrow(coords)} rows.
#'   When provided, \code{color_by} and \code{facet_by} can be column
#'   name strings, and \code{highlight} can be a named list of filter
#'   conditions applied to these columns.
#' @param facet_by Optional faceting variable(s). Either a character/factor
#'   vector for single-variable faceting (\code{facet_wrap}), or a
#'   data.frame with 1--2 columns for grid faceting (\code{facet_grid}).
#'   When \code{metadata} is provided, can be a character vector of 1--2
#'   column names.
#' @param highlight Optional. Which points to highlight, with
#'   non-highlighted points faded to a gray background. Accepts a logical
#'   vector, integer indices, or a named list for multi-column filtering
#'   (requires \code{metadata}; e.g.,
#'   \code{list(epitope = "PA", subject = c("S1", "S2"))}).
#' @param highlight_color Character. Color for highlighted points when
#'   \code{color_by} is \code{NULL}. Default \code{"#E41A1C"}.
#' @param background_alpha Numeric. Alpha for non-highlighted points.
#'   Default \code{0.15}.
#' @param tcr_rep A \code{\linkS4class{TCRrep}} object.  If provided and
#'   \code{metadata} is \code{NULL}, the clone data.frame is used as
#'   \code{metadata}, enabling column-name lookups for \code{color_by},
#'   \code{facet_by}, and \code{highlight}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' kpca <- compute_tcrdist_kernel_pca(tcr_df, "human")
#' plot_tcr_scatter(kpca$embeddings[, 1:2], color_by = tcr_df$epitope)
#'
#' # Faceting by epitope
#' plot_tcr_scatter(kpca$embeddings[, 1:2], color_by = tcr_df$epitope,
#'                  facet_by = tcr_df$epitope)
#'
#' # Highlight specific clones
#' plot_tcr_scatter(kpca$embeddings[, 1:2],
#'                  highlight = tcr_df$epitope == "PA")
#'
#' # Using metadata for column-name lookups
#' plot_tcr_scatter(kpca$embeddings[, 1:2],
#'                  metadata = tcr_df,
#'                  color_by = "epitope",
#'                  facet_by = c("epitope", "subject"),
#'                  highlight = list(epitope = "PA"))
#' }
#'
#' @seealso \code{\link{compute_tcrdist_kernel_pca}},
#'   \code{\link{compute_tcrdist_umap}},
#'   \code{\link{plot_tcrdist_heatmap}}
#' @export
plot_tcr_scatter <- function(coords, color_by = NULL, title = NULL,
                              point_size = 1, alpha = 1,
                              axis_label_prefix = "KPCA",
                              legend_title = NULL,
                              palette = NULL,
                              show_labels = FALSE,
                              label_size = 3,
                              na_color = "#DDDDDD",
                              metadata = NULL,
                              facet_by = NULL,
                              highlight = NULL,
                              highlight_color = "#E41A1C",
                              background_alpha = 0.15,
                              tcr_rep = NULL) {
    .check_ggplot2("plot_tcr_scatter()")

    if (!is.null(tcr_rep) && is.null(metadata)) {
        metadata <- .extract_from_tcr_rep(tcr_rep)$clone_df
    }

    coords <- as.matrix(coords)
    stopifnot(ncol(coords) == 2L)
    n <- nrow(coords)

    # --- Resolve metadata lookups early ---
    if (!is.null(metadata)) {
        stopifnot(is.data.frame(metadata), nrow(metadata) == n)

        if (!is.null(color_by) && is.character(color_by) &&
            length(color_by) == 1L && color_by %in% colnames(metadata)) {
            color_by <- metadata[[color_by]]
        }

        if (!is.null(facet_by) && is.character(facet_by) &&
            length(facet_by) <= 2L && all(facet_by %in% colnames(metadata))) {
            facet_by <- metadata[, facet_by, drop = FALSE]
        }
    }

    # --- Resolve highlight to logical vector ---
    if (!is.null(highlight)) {
        if (is.list(highlight) && !is.null(names(highlight))) {
            stopifnot(!is.null(metadata))
            hl <- rep(TRUE, n)
            for (nm in names(highlight)) {
                stopifnot(nm %in% colnames(metadata))
                hl <- hl & (metadata[[nm]] %in% highlight[[nm]])
            }
            highlight <- hl
        } else if (is.numeric(highlight) || is.integer(highlight)) {
            hl <- rep(FALSE, n)
            hl[highlight] <- TRUE
            highlight <- hl
        }
        stopifnot(is.logical(highlight), length(highlight) == n)
    }

    # --- Build internal data.frame ---
    df <- data.frame(
        x = coords[, 1L],
        y = coords[, 2L],
        stringsAsFactors = FALSE
    )

    # --- Attach facet columns ---
    facet_cols <- NULL
    if (!is.null(facet_by)) {
        if (is.data.frame(facet_by)) {
            stopifnot(nrow(facet_by) == n,
                      ncol(facet_by) >= 1L, ncol(facet_by) <= 2L)
            for (col in colnames(facet_by)) df[[col]] <- facet_by[[col]]
            facet_cols <- colnames(facet_by)
        } else {
            stopifnot(length(facet_by) == n)
            df$facet_var <- facet_by
            facet_cols <- "facet_var"
        }
    }

    # --- Build plot ---
    use_highlight <- !is.null(highlight)

    if (is.null(color_by)) {
        if (use_highlight) {
            bg_df <- df[!highlight, , drop = FALSE]
            fg_df <- df[highlight, , drop = FALSE]
            p <- ggplot2::ggplot(
                     mapping = ggplot2::aes(x = .data$x, y = .data$y)) +
                ggplot2::geom_point(data = bg_df, size = point_size,
                                    alpha = background_alpha,
                                    color = "gray80") +
                ggplot2::geom_point(data = fg_df, size = point_size,
                                    alpha = alpha, color = highlight_color)
        } else {
            p <- ggplot2::ggplot(df,
                     ggplot2::aes(x = .data$x, y = .data$y)) +
                ggplot2::geom_point(size = point_size, alpha = alpha,
                                    color = "gray50")
        }
    } else {
        stopifnot(length(color_by) == n)
        is_categorical <- is.factor(color_by) || is.character(color_by)
        df$color <- color_by

        if (use_highlight) {
            bg_df <- df[!highlight, , drop = FALSE]
            fg_df <- df[highlight, , drop = FALSE]
            if (!is_categorical) {
                fg_df <- fg_df[order(fg_df$color, na.last = FALSE), ]
            }
            p <- ggplot2::ggplot(
                     mapping = ggplot2::aes(x = .data$x, y = .data$y)) +
                ggplot2::geom_point(data = bg_df, size = point_size,
                                    alpha = background_alpha,
                                    color = "gray80") +
                ggplot2::geom_point(data = fg_df,
                                    ggplot2::aes(color = .data$color),
                                    size = point_size, alpha = alpha)
        } else {
            if (!is_categorical) {
                df <- df[order(df$color, na.last = FALSE), ]
            }
            p <- ggplot2::ggplot(df,
                     ggplot2::aes(x = .data$x, y = .data$y,
                                  color = .data$color)) +
                ggplot2::geom_point(size = point_size, alpha = alpha)
        }

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

            if (show_labels) {
                label_df <- if (use_highlight) fg_df else df
                if (!is.null(facet_cols)) {
                    grp_formula <- stats::as.formula(
                        paste("cbind(x, y) ~",
                              paste(c("color", facet_cols),
                                    collapse = " + ")))
                    centroids <- stats::aggregate(grp_formula,
                                                   data = label_df,
                                                   FUN = mean)
                } else {
                    centroids <- stats::aggregate(
                        cbind(x, y) ~ color, data = label_df, FUN = mean)
                }
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

    # --- Faceting ---
    if (!is.null(facet_cols)) {
        if (length(facet_cols) == 1L) {
            p <- p + ggplot2::facet_wrap(
                stats::as.formula(paste("~", facet_cols)))
        } else {
            p <- p + ggplot2::facet_grid(
                stats::as.formula(paste(facet_cols[1L], "~", facet_cols[2L])))
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
