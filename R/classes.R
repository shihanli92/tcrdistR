#' @include constants.R
NULL

#' TCRrep S4 class
#'
#' Represents a T-cell receptor repertoire with associated metadata, distance
#' matrices, and analysis results. This is the central data object of the
#' tcrdistR package. Create instances using the \code{\link{TCRrep}} constructor
#' rather than calling \code{new("TCRrep", ...)} directly.
#'
#' @slot clone_df A \code{data.frame} of clonotypes. Required columns depend on
#'   the \code{chains} value: alpha-beta (\code{"AB"}) requires \code{va},
#'   \code{cdr3a}, \code{vb}, \code{cdr3b}; alpha only (\code{"A"}) requires
#'   \code{va}, \code{cdr3a}; beta only (\code{"B"}) requires \code{vb},
#'   \code{cdr3b}; gamma-delta (\code{"GD"}) requires \code{va}, \code{cdr3a},
#'   \code{vb}, \code{cdr3b}.
#' @slot organism Character string specifying the organism. Must be one of the
#'   organisms available in the bundled gene database (e.g. \code{"human"},
#'   \code{"mouse"}).
#' @slot chains Character string specifying the chain combination. One of
#'   \code{"AB"}, \code{"A"}, \code{"B"}, or \code{"GD"}.
#' @slot metric Character string specifying the distance metric. One of
#'   \code{"tcrdist"}, \code{"hamming"}, \code{"levenshtein"}, or \code{"nw"}.
#' @slot weights A named list with elements \code{cdr3} and \code{v_region}
#'   specifying the integer weights applied to each region's distance
#'   contribution.
#' @slot gap_penalties A named list with elements \code{cdr3} and \code{v_region}
#'   specifying the integer gap penalties used in sequence alignment.
#' @slot dist_a Distance matrix for alpha-chain only comparisons. Either
#'   \code{NULL} or a numeric matrix.
#' @slot dist_b Distance matrix for beta-chain only comparisons. Either
#'   \code{NULL} or a numeric matrix.
#' @slot paired_dist Paired (alpha + beta combined) distance matrix. Either
#'   \code{NULL} or a numeric matrix. Populated when
#'   \code{compute_distances = TRUE} is passed to \code{\link{TCRrep}}.
#' @slot knn_indices K-nearest-neighbour index matrix. Either \code{NULL} or
#'   an integer matrix of dimensions N x K.
#' @slot knn_distances K-nearest-neighbour distance matrix. Either \code{NULL}
#'   or a numeric matrix of dimensions N x K.
#' @slot rep_dists A named list of representative-clonotype distance objects,
#'   populated by downstream analyses.
#' @slot meta_clonotypes Meta-clonotype table. Either \code{NULL} or a
#'   \code{data.frame}.
#' @slot motifs A named list of motif objects, populated by downstream analyses.
#'
#' @seealso \code{\link{TCRrep}} for the recommended constructor.
#'
#' @import methods
#' @exportClass TCRrep
setClass(
    "TCRrep",
    representation(
        clone_df      = "data.frame",
        organism      = "character",
        chains        = "character",
        metric        = "character",
        weights       = "list",
        gap_penalties = "list",
        dist_a        = "ANY",
        dist_b        = "ANY",
        paired_dist   = "ANY",
        knn_indices   = "ANY",
        knn_distances = "ANY",
        rep_dists     = "list",
        meta_clonotypes = "ANY",
        motifs        = "list"
    ),
    prototype(
        clone_df        = data.frame(),
        organism        = "human",
        chains          = "AB",
        metric          = "tcrdist",
        weights         = list(cdr3 = 3L, v_region = 1L),
        gap_penalties   = list(cdr3 = 12L, v_region = 4L),
        dist_a          = NULL,
        dist_b          = NULL,
        paired_dist     = NULL,
        knn_indices     = NULL,
        knn_distances   = NULL,
        rep_dists       = list(),
        meta_clonotypes = NULL,
        motifs          = list()
    )
)

# ---------------------------------------------------------------------------
# Validity
# ---------------------------------------------------------------------------

setValidity("TCRrep", function(object) {
    errors <- character(0L)

    # ---- organism ------------------------------------------------------------
    if (length(object@organism) != 1L || !nzchar(object@organism)) {
        errors <- c(errors,
            "'organism' must be a non-empty character string of length 1")
    }

    # ---- chains --------------------------------------------------------------
    valid_chains <- c("AB", "A", "B", "GD")
    if (length(object@chains) != 1L ||
        !(object@chains %in% valid_chains)) {
        errors <- c(errors, sprintf(
            "'chains' must be one of: %s",
            paste(valid_chains, collapse = ", ")
        ))
    }

    # ---- metric --------------------------------------------------------------
    valid_metrics <- c("tcrdist", "hamming", "levenshtein", "nw")
    if (length(object@metric) != 1L ||
        !(object@metric %in% valid_metrics)) {
        errors <- c(errors, sprintf(
            "'metric' must be one of: %s",
            paste(valid_metrics, collapse = ", ")
        ))
    }

    # ---- clone_df column requirements (only checked when non-empty) ---------
    if (nrow(object@clone_df) > 0L) {
        chains_val <- object@chains
        required_cols <- switch(
            chains_val,
            "AB" = c("va", "cdr3a", "vb", "cdr3b"),
            "A"  = c("va", "cdr3a"),
            "B"  = c("vb", "cdr3b"),
            "GD" = c("va", "cdr3a", "vb", "cdr3b"),
            character(0L)
        )
        missing_cols <- setdiff(required_cols, colnames(object@clone_df))
        if (length(missing_cols) > 0L) {
            errors <- c(errors, sprintf(
                "'clone_df' is missing required columns for chains='%s': %s",
                chains_val,
                paste(missing_cols, collapse = ", ")
            ))
        }
    }

    # ---- weights -------------------------------------------------------------
    required_weight_names <- c("cdr3", "v_region")
    missing_weights <- setdiff(required_weight_names, names(object@weights))
    if (length(missing_weights) > 0L) {
        errors <- c(errors, sprintf(
            "'weights' must contain named elements: %s (missing: %s)",
            paste(required_weight_names, collapse = ", "),
            paste(missing_weights, collapse = ", ")
        ))
    }

    # ---- gap_penalties -------------------------------------------------------
    required_gap_names <- c("cdr3", "v_region")
    missing_gaps <- setdiff(required_gap_names, names(object@gap_penalties))
    if (length(missing_gaps) > 0L) {
        errors <- c(errors, sprintf(
            "'gap_penalties' must contain named elements: %s (missing: %s)",
            paste(required_gap_names, collapse = ", "),
            paste(missing_gaps, collapse = ", ")
        ))
    }

    if (length(errors) == 0L) TRUE else errors
})
