#' @include classes.R
NULL

#' Create a TCRrep object
#'
#' Constructs a \code{\link{TCRrep}} S4 object from a clonotype data frame and
#' optional parameters. This is the recommended way to create a \code{TCRrep}
#' instance (following Bioconductor convention, rather than calling
#' \code{new("TCRrep", ...)} directly).
#'
#' Factor columns (\code{va}, \code{cdr3a}, \code{vb}, \code{cdr3b}) are
#' automatically coerced to character. The organism is validated against the
#' bundled gene database on construction.
#'
#' @param clone_df A \code{data.frame} of clonotypes. Required columns depend
#'   on the \code{chains} argument:
#'   \describe{
#'     \item{\code{"AB"}}{Requires \code{va}, \code{cdr3a}, \code{vb},
#'       \code{cdr3b}.}
#'     \item{\code{"A"}}{Requires \code{va}, \code{cdr3a}.}
#'     \item{\code{"B"}}{Requires \code{vb}, \code{cdr3b}.}
#'     \item{\code{"GD"}}{Requires \code{va}, \code{cdr3a}, \code{vb},
#'       \code{cdr3b}.}
#'   }
#' @param organism Character string. Organism key recognised by
#'   \code{\link{load_gene_database}}, e.g. \code{"human"} or \code{"mouse"}.
#' @param chains Character string. One of \code{"AB"} (default), \code{"A"},
#'   \code{"B"}, or \code{"GD"}.
#' @param metric Character string. Distance metric to use. One of
#'   \code{"tcrdist"} (default), \code{"hamming"}, \code{"levenshtein"}, or
#'   \code{"nw"}.
#' @param compute_distances Logical. If \code{TRUE} and \code{nrow(clone_df) > 0},
#'   compute the pairwise distance matrix immediately and store it in the
#'   \code{paired_dist} slot. Defaults to \code{FALSE}.
#' @param weight_cdr3 Integer. Weight applied to CDR3 distances. Defaults to
#'   \code{WEIGHT_CDR3_REGION} (3L).
#' @param gap_penalty_cdr3 Integer. Gap penalty for CDR3 alignments. Defaults
#'   to \code{GAP_PENALTY_CDR3_REGION} (12L).
#' @param weight_v_region Integer. Weight applied to V-region distances.
#'   Defaults to \code{WEIGHT_V_REGION} (1L).
#' @param gap_penalty_v_region Integer. Gap penalty for V-region alignments.
#'   Defaults to \code{GAP_PENALTY_V_REGION} (4L).
#'
#' @return A valid \code{\link{TCRrep}} S4 object.
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
#'
#' # Basic construction (no distance computation)
#' obj <- TCRrep(tcrs, organism = "human")
#'
#' # With distance computation
#' obj <- TCRrep(tcrs, organism = "human", compute_distances = TRUE)
#' dim(obj@paired_dist)  # 2 x 2
#'
#' # Alpha chain only
#' obj_a <- TCRrep(tcrs[, c("va", "cdr3a")], organism = "human", chains = "A")
#' }
#'
#' @export
TCRrep <- function(clone_df,
                   organism             = "human",
                   chains               = "AB",
                   metric               = "tcrdist",
                   compute_distances    = FALSE,
                   weight_cdr3          = WEIGHT_CDR3_REGION,
                   gap_penalty_cdr3     = GAP_PENALTY_CDR3_REGION,
                   weight_v_region      = WEIGHT_V_REGION,
                   gap_penalty_v_region = GAP_PENALTY_V_REGION) {

    # ---- Input validation ---------------------------------------------------
    if (!is.data.frame(clone_df)) {
        stop("TCRrep: 'clone_df' must be a data.frame")
    }

    # ---- Coerce factor columns to character ---------------------------------
    factor_cols <- c("va", "cdr3a", "vb", "cdr3b")
    for (col in factor_cols) {
        if (col %in% colnames(clone_df) && is.factor(clone_df[[col]])) {
            clone_df[[col]] <- as.character(clone_df[[col]])
        }
    }

    # ---- Validate organism against the gene database ------------------------
    # This call errors with an informative message if organism is unknown.
    load_gene_database(organism)

    # ---- Assemble weights and gap_penalties lists ----------------------------
    weights_list <- list(
        cdr3     = as.integer(weight_cdr3),
        v_region = as.integer(weight_v_region)
    )
    gap_penalties_list <- list(
        cdr3     = as.integer(gap_penalty_cdr3),
        v_region = as.integer(gap_penalty_v_region)
    )

    # ---- Construct the S4 object (triggers validity check) ------------------
    obj <- new(
        "TCRrep",
        clone_df      = clone_df,
        organism      = organism,
        chains        = chains,
        metric        = metric,
        weights       = weights_list,
        gap_penalties = gap_penalties_list
    )

    # ---- Optionally compute distances ---------------------------------------
    if (isTRUE(compute_distances) && nrow(clone_df) > 0L) {
        obj@paired_dist <- tcrdist_matrix(
            clone_df,
            organism,
            weight_cdr3      = weight_cdr3,
            gap_penalty_cdr3 = gap_penalty_cdr3
        )
    }

    obj
}
