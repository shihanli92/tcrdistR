# Shared plotting utilities for tcrdistR visualization functions
#
# Internal helpers: palette generation, dependency checks.
# No exported functions in this file.


# Suppress R CMD check NOTEs for ggplot2 .data pronoun
utils::globalVariables(".data")


# ---------------------------------------------------------------------------
# .check_ggplot2
# ---------------------------------------------------------------------------

#' Check that ggplot2 is available
#'
#' @param fn Character. Name of the calling function (for the error message).
#' @return Invisible \code{TRUE} if available; otherwise stops with an
#'   informative error.
#' @keywords internal
#' @noRd
.check_ggplot2 <- function(fn = "this function") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop(sprintf(
            "Package 'ggplot2' is required for %s. Install it with: install.packages(\"ggplot2\")",
            fn
        ), call. = FALSE)
    }
    invisible(TRUE)
}


# ---------------------------------------------------------------------------
# .tcrdistR_palette
# ---------------------------------------------------------------------------

#' Generate a qualitative color palette matching matplotlib tab10/tab20
#'
#' Returns \code{n} colors from the matplotlib tab10 (for n <= 10) or tab20
#' (for n > 10) palettes, recycled as needed.
#'
#' @param n Integer. Number of colors needed.
#' @return A character vector of \code{n} hex color strings.
#' @keywords internal
#' @noRd
.tcrdistR_palette <- function(n) {
    tab10 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
    tab20 <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
               "#2ca02c", "#98df8a", "#d62728", "#ff9896",
               "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
               "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7",
               "#bcbd22", "#dbdb8d", "#17becf", "#9edae5")
    if (n <= 10L) {
        rep_len(tab10, n)
    } else {
        rep_len(tab20, n)
    }
}
