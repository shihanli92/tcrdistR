# Distance-based fuzzy joining of two TCR datasets
#
# Ported from Python tcrdist3 join.py.


# ---------------------------------------------------------------------------
# tcrdist_join  (exported)
# ---------------------------------------------------------------------------

#' Join two TCR datasets by TCRdist threshold
#'
#' Performs a fuzzy merge of two TCR data.frames: for each row in
#' \code{left_df}, finds rows in \code{right_df} within \code{radius}
#' TCRdist units. Similar to a SQL JOIN but using distance instead of
#' exact matching.
#'
#' @param left_df Data.frame with TCR columns (\code{va}, \code{vb},
#'   \code{cdr3a}, \code{cdr3b}) plus any extra columns.
#' @param right_df Data.frame with TCR columns (\code{va}, \code{vb},
#'   \code{cdr3a}, \code{cdr3b}) plus any extra columns.
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#' @param radius Numeric. Maximum TCRdist for a match.
#' @param max_n Integer. Maximum number of matches per left row. Default
#'   \code{5L}. Closest matches are kept.
#' @param type Character string. Join type: \code{"inner"} (default) keeps
#'   only matched pairs, \code{"left"} keeps all left rows (NAs for
#'   unmatched).
#' @param suffix Character vector of length 2. Suffixes for disambiguating
#'   column names. Default \code{c("_x", "_y")}.
#'
#' @return A data.frame with columns from both sides (suffixed if
#'   overlapping) plus a \code{tcrdist} column.
#'
#' @examples
#' \dontrun{
#' matches <- tcrdist_join(query_tcrs, reference_tcrs, "human", radius = 50)
#' }
#'
#' @export
tcrdist_join <- function(left_df, right_df, organism, radius,
                          max_n = 5L, type = c("inner", "left"),
                          suffix = c("_x", "_y")) {
    type <- match.arg(type)
    stopifnot(
        is.data.frame(left_df), is.data.frame(right_df),
        is.numeric(radius), length(radius) == 1L, radius >= 0
    )

    n_left <- nrow(left_df)
    n_right <- nrow(right_df)

    if (n_left == 0L || n_right == 0L) {
        return(.empty_join_result(left_df, right_df, suffix))
    }

    # Compute rectangular distances (dense, with threshold applied)
    rect_mat <- tcrdist_rect(left_df, right_df, organism)

    # Build result rows
    left_names <- paste0(colnames(left_df), suffix[1L])
    right_names <- paste0(colnames(right_df), suffix[2L])
    result_rows <- vector("list", n_left)  # pre-allocate

    for (i in seq_len(n_left)) {
        dists <- rect_mat[i, ]
        within <- which(dists <= radius)

        if (length(within) == 0L) {
            if (type == "left") {
                left_row <- left_df[i, , drop = FALSE]
                colnames(left_row) <- left_names
                right_row <- right_df[1L, , drop = FALSE]
                right_row[1L, ] <- NA
                colnames(right_row) <- right_names
                result_rows[[i]] <- cbind(left_row, right_row,
                                          data.frame(tcrdist = NA_real_),
                                          row.names = NULL)
            }
            next
        }

        # Sort by distance, keep top max_n
        ord <- order(dists[within])
        if (length(ord) > max_n) ord <- ord[seq_len(max_n)]
        sel <- within[ord]

        left_row <- left_df[rep(i, length(sel)), , drop = FALSE]
        colnames(left_row) <- left_names
        right_row <- right_df[sel, , drop = FALSE]
        colnames(right_row) <- right_names
        result_rows[[i]] <- cbind(left_row, right_row,
                                  data.frame(tcrdist = dists[sel]),
                                  row.names = NULL)
    }

    result_rows <- result_rows[!vapply(result_rows, is.null, logical(1L))]

    if (length(result_rows) == 0L) {
        return(.empty_join_result(left_df, right_df, suffix))
    }

    result <- do.call(rbind, result_rows)
    rownames(result) <- NULL
    result
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.empty_join_result <- function(left_df, right_df, suffix) {
    left_cols <- paste0(colnames(left_df), suffix[1L])
    right_cols <- paste0(colnames(right_df), suffix[2L])
    all_cols <- c(left_cols, right_cols, "tcrdist")

    empty <- data.frame(matrix(nrow = 0, ncol = length(all_cols)))
    colnames(empty) <- all_cols
    empty
}
