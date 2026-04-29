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
#' @param tcr_df Data.frame with TCR columns.
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#' @param method Clustering method for \code{stats::hclust()}. Default
#'   \code{"average"} (UPGMA).
#' @param max_tcrs Integer. Subsample if N exceeds this. Default \code{2000L}.
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
tcrdist_hclust <- function(tcr_df, organism, method = "average",
                            max_tcrs = 2000L) {
    n <- nrow(tcr_df)
    stopifnot(n >= 2L)

    indices <- seq_len(n)
    if (n > max_tcrs) {
        indices <- sort(sample(n, max_tcrs))
        tcr_df <- tcr_df[indices, , drop = FALSE]
    }

    dist_mat <- tcrdist_matrix(tcr_df, organism)
    hc <- stats::hclust(stats::as.dist(dist_mat), method = method)

    list(
        hclust = hc,
        dist_matrix = dist_mat,
        indices = indices
    )
}


# ---------------------------------------------------------------------------
# cluster_tcrs  (exported)
# ---------------------------------------------------------------------------

#' Cluster TCRs by TCRdist-based hierarchical clustering
#'
#' Computes TCRdist distances, performs hierarchical clustering, and cuts
#' the tree into groups using either a fixed number of clusters (\code{k})
#' or a height threshold (\code{h}).
#'
#' @param tcr_df Data.frame with TCR columns.
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#' @param k Integer. Number of clusters. Exactly one of \code{k} or
#'   \code{h} must be specified.
#' @param h Numeric. Height at which to cut the dendrogram. Exactly one
#'   of \code{k} or \code{h} must be specified.
#' @param method Clustering method. Default \code{"average"}.
#'
#' @return An integer vector of cluster assignments (length
#'   \code{nrow(tcr_df)}).
#'
#' @examples
#' \dontrun{
#' clusters <- cluster_tcrs(tcr_df, "human", k = 5)
#' table(clusters)
#' }
#'
#' @seealso \code{\link{tcrdist_hclust}}, \code{\link{neighborhood_test}}
#' @export
cluster_tcrs <- function(tcr_df, organism, k = NULL, h = NULL,
                          method = "average") {
    if (is.null(k) && is.null(h)) {
        stop("cluster_tcrs: exactly one of 'k' or 'h' must be specified",
             call. = FALSE)
    }
    if (!is.null(k) && !is.null(h)) {
        stop("cluster_tcrs: specify either 'k' or 'h', not both",
             call. = FALSE)
    }

    result <- tcrdist_hclust(tcr_df, organism, method = method)
    stats::cutree(result$hclust, k = k, h = h)
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
#' @param tcr_df Data.frame with TCR columns.
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#' @param variable Character or factor vector of length \code{nrow(tcr_df)}.
#'   The categorical variable to test.
#' @param radius Numeric. Maximum TCRdist for neighborhood membership.
#'   Default \code{50}.
#' @param test Character string. \code{"fisher"} (default, for binary) or
#'   \code{"chisq"} (for multi-category).
#' @param p_adjust_method Character string. Method for \code{stats::p.adjust()}.
#'   Default \code{"BH"} (Benjamini-Hochberg).
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
neighborhood_test <- function(tcr_df, organism, variable, radius = 50,
                                test = c("fisher", "chisq"),
                                p_adjust_method = "BH") {
    test <- match.arg(test)
    n <- nrow(tcr_df)
    stopifnot(length(variable) == n)

    variable <- as.factor(variable)
    lvls <- levels(variable)
    n_lvls <- length(lvls)

    if (test == "fisher" && n_lvls != 2L) {
        stop("neighborhood_test: Fisher's exact test requires exactly 2 ",
             "categories. Use test='chisq' for ", n_lvls, " categories.",
             call. = FALSE)
    }

    # Compute dense distance matrix
    dist_mat <- tcrdist_matrix(tcr_df, organism)

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
