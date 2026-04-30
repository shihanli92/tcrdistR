#' @include classes.R
NULL

# ---------------------------------------------------------------------------
# Subsetting: [
# ---------------------------------------------------------------------------

#' Subset a TCRrep object by row indices
#'
#' Subsets a \code{\link{TCRrep}} by clonotype indices. Filters
#' \code{clone_df} and all associated distance matrices simultaneously.
#' KNN matrices and meta-clonotypes are cleared since they reference the
#' original indexing.
#'
#' @param x A \code{TCRrep} object.
#' @param i Row indices: logical, integer, or character (rownames).
#' @param j Ignored.
#' @param ... Ignored.
#' @param drop Ignored.
#' @return A new \code{TCRrep} object with the selected clonotypes.
#'
#' @examples
#' \donttest{
#' data(dash)
#' rep <- TCRrep(dash, organism = "mouse", compute_distances = TRUE)
#' pa <- rep[rep@@clone_df$epitope == "PA", ]
#' pa
#' }
#'
#' @name subset-TCRrep
#' @aliases [,TCRrep,ANY,ANY,ANY-method
#' @export
setMethod("[", signature(x = "TCRrep"), function(x, i, j, ..., drop = TRUE) {
    if (is.logical(i)) {
        if (length(i) != nrow(x@clone_df)) {
            stop(sprintf(
                "logical index length (%d) must match number of clonotypes (%d)",
                length(i), nrow(x@clone_df)), call. = FALSE)
        }
        i <- which(i)
    }

    obj <- x
    obj@clone_df <- x@clone_df[i, , drop = FALSE]
    if (!is.null(x@paired_dist)) obj@paired_dist <- x@paired_dist[i, i, drop = FALSE]
    if (!is.null(x@dist_a))      obj@dist_a      <- x@dist_a[i, i, drop = FALSE]
    if (!is.null(x@dist_b))      obj@dist_b      <- x@dist_b[i, i, drop = FALSE]
    obj@knn_indices     <- NULL
    obj@knn_distances   <- NULL
    obj@meta_clonotypes <- NULL
    obj
})


# ---------------------------------------------------------------------------
# subset
# ---------------------------------------------------------------------------

#' Subset a TCRrep using non-standard evaluation
#'
#' Filter clonotypes using column expressions evaluated against
#' \code{clone_df}. This is the recommended way to filter a TCRrep:
#'
#' \preformatted{rep |> subset(epitope == "PA")}
#'
#' @param x A \code{TCRrep} object.
#' @param subset A logical expression evaluated in the context of
#'   \code{x@@clone_df}. Column names can be used directly.
#' @param ... Ignored.
#' @return A new \code{TCRrep} object with only matching clonotypes.
#'
#' @examples
#' \donttest{
#' data(dash)
#' rep <- TCRrep(dash, organism = "mouse", compute_distances = TRUE)
#' pa <- subset(rep, epitope == "PA")
#' pa
#' }
#'
#' @export
setMethod("subset", "TCRrep", function(x, subset, ...) {
    expr <- substitute(subset)
    # S4 dispatch adds frames; parent.frame(2) reaches the actual caller
    idx <- eval(expr, x@clone_df, parent.frame(2))
    if (!is.logical(idx)) {
        stop("'subset' must evaluate to a logical vector", call. = FALSE)
    }
    idx[is.na(idx)] <- FALSE
    x[idx, ]
})


# ---------------------------------------------------------------------------
# show
# ---------------------------------------------------------------------------

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
