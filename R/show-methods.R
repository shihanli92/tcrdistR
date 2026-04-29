#' @include classes.R
NULL

#' Show method for TCRrep objects
#'
#' Prints a compact summary of a \code{\link{TCRrep}} object, including the
#' number of clonotypes, organism, chain combination, distance metric, distance
#' computation status, and k-nearest-neighbour status when available.
#'
#' @param object A \code{TCRrep} object.
#' @return Returns \code{invisible(object)}.
#'
#' @examples
#' \donttest{
#' tcrs <- data.frame(
#'     va    = c("TRAV1-1*01", "TRAV1-1*01"),
#'     cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
#'     vb    = c("TRBV19*01", "TRBV19*01"),
#'     cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
#'     stringsAsFactors = FALSE
#' )
#' obj <- TCRrep(tcrs, organism = "human")
#' show(obj)
#' }
#'
#' @export
setMethod("show", "TCRrep", function(object) {
    n_clonotypes <- nrow(object@clone_df)
    cat(sprintf("TCRrep object: %d clonotypes\n", n_clonotypes))
    cat(sprintf("  organism: %s\n", object@organism))
    cat(sprintf("  chains: %s\n",   object@chains))
    cat(sprintf("  metric: %s\n",   object@metric))

    # ---- Distance status ----------------------------------------------------
    dist_status <- if (!is.null(object@paired_dist)) {
        dist_type <- if (inherits(object@paired_dist, "sparseMatrix")) {
            "sparse"
        } else {
            "dense"
        }
        sprintf("paired(%s)", dist_type)
    } else if (!is.null(object@dist_a) || !is.null(object@dist_b)) {
        parts <- character(0L)
        if (!is.null(object@dist_a)) parts <- c(parts, "alpha")
        if (!is.null(object@dist_b)) parts <- c(parts, "beta")
        paste(parts, collapse = "+")
    } else {
        "not computed"
    }
    cat(sprintf("  distances: %s\n", dist_status))

    # ---- KNN status ---------------------------------------------------------
    if (!is.null(object@knn_indices)) {
        k_val <- ncol(object@knn_indices)
        cat(sprintf("  knn: K=%d\n", k_val))
    }

    invisible(object)
})
