#' tcrdistR: C++-Accelerated TCR Distance Calculations
#'
#' Compute pairwise TCRdist distances for T-cell receptor (TCR) repertoire
#' analysis. All distance computations are implemented in C++ via Rcpp for
#' high performance, with a thin R layer for input validation, I/O, and
#' visualization. Targets feature parity with Python tcrdist3.
#'
#' @details
#' The package provides the following functional areas:
#'
#' \strong{Distance computation:}
#' \code{\link{tcrdist_matrix}}, \code{\link{tcrdist_sparse}},
#' \code{\link{tcrdist_rect}}, \code{\link{weighted_cdr3_distance}},
#' \code{\link{bsd4_matrix}}, \code{\link{hamming_distance}},
#' \code{\link{hamming_matrix}}
#'
#' \strong{Neighbor search:}
#' \code{\link{tcrdist_knn}}, \code{\link{tcrdist_radius_neighbors}},
#' \code{\link{knn_from_matrix}}, \code{\link{knn_from_pca}}
#'
#' \strong{Clumping and background models:}
#' \code{\link{find_clumping}}, \code{\link{setup_tcr_groups}}
#'
#' \strong{Database matching:}
#' \code{\link{find_significant_tcrdist_matches}},
#' \code{\link{match_tcrs_to_db}},
#' \code{\link{strict_single_chain_match_tcrs_to_db}}
#'
#' \strong{Dimensionality reduction:}
#' \code{\link{compute_tcrdist_kernel_pca}}
#'
#' \strong{CD8 scoring:}
#' \code{\link{make_cd8_score_table_column}}
#'
#' \strong{Diversity metrics:}
#' \code{\link{tcr_diversity}}, \code{\link{tcr_fuzzy_diversity}},
#' \code{\link{tcr_richness}}, \code{\link{tcr_clonality}}
#'
#' \strong{Clustering and neighborhood tests:}
#' \code{\link{tcrdist_hclust}}, \code{\link{cluster_tcrs}},
#' \code{\link{neighborhood_test}}
#'
#' \strong{Meta-clonotypes:}
#' \code{\link{find_meta_clonotypes}},
#' \code{\link{summarize_meta_clonotype}}
#'
#' \strong{Joins:}
#' \code{\link{tcrdist_join}}
#'
#' \strong{Visualization:}
#' \code{\link{plot_tcrdist_heatmap}}, \code{\link{plot_tcrdist_dendrogram}},
#' \code{\link{plot_distance_distribution}}, \code{\link{plot_cdr3_logo}},
#' \code{\link{plot_junction_bars}}, \code{\link{plot_gene_usage}},
#' \code{\link{plot_tcr_scatter}}
#'
#' \strong{I/O:}
#' \code{\link{read_tcr_table}}, \code{\link{read_airr}},
#' \code{\link{read_adaptive}}, \code{\link{read_10x}}
#'
#' @seealso \code{vignette("tcrdistR-getting-started")} for a tutorial
#'   introduction.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib tcrdistR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import methods
#' @importFrom Matrix sparseMatrix
NULL

# Package-private environment for caching
.tcrdistR_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
    .tcrdistR_env$all_genes <- list()
    .tcrdistR_env$rep_dists <- list()
    .tcrdistR_env$v_dist_matrices <- list()
}
