#' @include classes.R
NULL

# ---------------------------------------------------------------------------
# .deduplicate_clones (internal helper)
# ---------------------------------------------------------------------------

#' Deduplicate clones by grouping columns
#'
#' Groups rows of \code{clone_df} by \code{group_cols}, sums the \code{count}
#' column for merged rows (adding \code{count = 1L} if absent), drops the
#' \code{clone_id} column (no longer meaningful after merge), and removes
#' rows with \code{NA} in any grouping column.
#'
#' @param clone_df A data.frame of clonotypes.
#' @param group_cols Character vector of column names to group by.
#' @return A deduplicated data.frame.
#' @noRd
.deduplicate_clones <- function(clone_df, group_cols) {
    missing <- setdiff(group_cols, colnames(clone_df))
    if (length(missing) > 0L) {
        stop(sprintf("deduplicate: columns not found in clone_df: %s",
                     paste(missing, collapse = ", ")), call. = FALSE)
    }

    # Drop rows with NA in any grouping column
    complete <- complete.cases(clone_df[, group_cols, drop = FALSE])
    if (!all(complete)) {
        n_dropped <- sum(!complete)
        message(sprintf(
            "deduplicate: dropping %d row(s) with NA in grouping columns",
            n_dropped
        ))
        clone_df <- clone_df[complete, , drop = FALSE]
    }

    if (nrow(clone_df) == 0L) return(clone_df)

    # Build composite key for grouping
    dup_key <- do.call(paste, c(clone_df[, group_cols, drop = FALSE],
                                list(sep = "\x1f")))
    if (!anyDuplicated(dup_key)) return(clone_df)

    # Add count column if missing
    has_count <- "count" %in% colnames(clone_df)
    if (!has_count) clone_df$count <- 1L

    # Aggregate: keep first row per group, sum counts
    n_before <- nrow(clone_df)
    split_idx <- split(seq_len(n_before), dup_key)
    result_list <- lapply(split_idx, function(idx) {
        row <- clone_df[idx[1L], , drop = FALSE]
        if (length(idx) > 1L) {
            row$count <- sum(clone_df$count[idx])
        }
        row
    })
    clone_df <- do.call(rbind, result_list)
    rownames(clone_df) <- NULL

    # Drop clone_id (no longer meaningful after merge)
    if ("clone_id" %in% colnames(clone_df)) {
        clone_df$clone_id <- NULL
    }

    message(sprintf("deduplicate: %d -> %d clones", n_before, nrow(clone_df)))
    clone_df
}


# ---------------------------------------------------------------------------
# TCRrep constructor (exported)
# ---------------------------------------------------------------------------

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
#' @param deduplicate Controls clone deduplication (matching tcrdist3 behavior).
#'   \describe{
#'     \item{\code{TRUE} (default)}{Deduplicate using chain columns
#'       (\code{va}, \code{cdr3a}, \code{vb}, \code{cdr3b}) plus \code{subject}
#'       if present. Within-subject duplicates are merged and \code{count}
#'       values summed.}
#'     \item{\code{FALSE}}{No deduplication; \code{clone_df} is stored as-is.}
#'     \item{Character vector}{Custom grouping columns. Only rows identical
#'       across all specified columns are merged. Example:
#'       \code{c("va", "cdr3a", "vb", "cdr3b")} to ignore subject.}
#'   }
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
#' # Basic construction (deduplicates by default)
#' obj <- TCRrep(tcrs, organism = "human")
#'
#' # With distance computation
#' obj <- TCRrep(tcrs, organism = "human", compute_distances = TRUE)
#' dim(obj@paired_dist)  # 2 x 2
#'
#' # Custom dedup columns (ignore subject, collapse across individuals)
#' obj <- TCRrep(tcrs, organism = "human",
#'               deduplicate = c("va", "cdr3a", "vb", "cdr3b"))
#'
#' # No deduplication
#' obj <- TCRrep(tcrs, organism = "human", deduplicate = FALSE)
#' }
#'
#' @seealso \code{\link{tcrdist_matrix}}, \code{\link{read_tcr_table}}
#' @export
TCRrep <- function(clone_df,
                   organism             = "human",
                   chains               = "AB",
                   deduplicate          = TRUE,
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

    # ---- Deduplication ------------------------------------------------------
    if (!isFALSE(deduplicate) && nrow(clone_df) > 0L) {
        if (is.character(deduplicate)) {
            group_cols <- deduplicate
        } else {
            chain_cols <- switch(chains,
                "AB" = c("va", "cdr3a", "vb", "cdr3b"),
                "A"  = c("va", "cdr3a"),
                "B"  = c("vb", "cdr3b"),
                "GD" = c("va", "cdr3a", "vb", "cdr3b"),
                character(0L)
            )
            group_cols <- chain_cols
            if ("subject" %in% colnames(clone_df)) {
                group_cols <- c(group_cols, "subject")
            }
        }
        clone_df <- .deduplicate_clones(clone_df, group_cols)
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
