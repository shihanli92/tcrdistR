#' @keywords internal
"_PACKAGE"

#' @useDynLib tcrdistR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Package-private environment for caching
.tcrdistR_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
    .tcrdistR_env$all_genes <- list()
    .tcrdistR_env$rep_dists <- list()
    .tcrdistR_env$v_dist_matrices <- list()
}
