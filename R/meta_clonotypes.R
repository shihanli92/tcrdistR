# Meta-clonotype detection: identify quasi-public TCR motifs
#
# Ported from Python tcrdist3 public.py and centers.py.


# ---------------------------------------------------------------------------
# .compute_distance_ecdf  (internal)
# ---------------------------------------------------------------------------

#' Compute empirical CDF of distances
#'
#' @param dists Numeric vector. Distances to background TCRs.
#' @param weights Numeric vector. Weights for each background TCR (optional).
#' @return A function mapping distance threshold -> proportion of background
#'   within that threshold.
#' @keywords internal
#' @noRd
.compute_distance_ecdf <- function(dists, weights = NULL) {
    dists <- dists[!is.na(dists)]
    if (length(dists) == 0L) {
        return(function(d) rep(0, length(d)))
    }

    if (is.null(weights)) {
        stats::ecdf(dists)
    } else {
        weights <- weights[!is.na(dists)]
        total_w <- sum(weights)
        if (total_w == 0) return(function(d) rep(0, length(d)))

        ord <- order(dists)
        sorted_d <- dists[ord]
        sorted_w <- weights[ord]
        cum_w <- cumsum(sorted_w) / total_w

        function(d) {
            vapply(d, function(x) {
                idx <- findInterval(x, sorted_d)
                if (idx == 0L) 0 else cum_w[idx]
            }, numeric(1L))
        }
    }
}


# ---------------------------------------------------------------------------
# .calc_radii  (internal)
# ---------------------------------------------------------------------------

#' Calculate per-TCR optimal radii from background ECDF
#'
#' For each enriched TCR, finds the maximum distance threshold at which
#' the proportion of background TCRs within that distance does not exceed
#' \code{ctrl_bkgd}.
#'
#' @param tcr_df Data.frame. Enriched TCRs.
#' @param organism Character string.
#' @param background_df Data.frame. Background TCRs.
#' @param ctrl_bkgd Numeric. Maximum background proportion allowed.
#' @param max_radius Numeric. Maximum radius to consider.
#' @return Integer vector of per-TCR radii.
#' @keywords internal
#' @noRd
.calc_radii <- function(tcr_df, organism, background_df,
                         ctrl_bkgd = 1e-5, max_radius = 50L) {
    n_enriched <- nrow(tcr_df)

    # Compute rectangular distances: enriched (rows) vs background (cols)
    rect_mat <- tcrdist_rect(tcr_df, background_df, organism)

    radii <- integer(n_enriched)
    test_thresholds <- seq(0L, max_radius, by = 1L)

    for (i in seq_len(n_enriched)) {
        ecdf_fn <- .compute_distance_ecdf(rect_mat[i, ])
        ecdf_vals <- ecdf_fn(test_thresholds)

        # Find maximum threshold where ECDF <= ctrl_bkgd
        valid <- which(ecdf_vals <= ctrl_bkgd)
        if (length(valid) == 0L) {
            radii[i] <- 0L
        } else {
            radii[i] <- test_thresholds[max(valid)]
        }
    }

    radii
}


# ---------------------------------------------------------------------------
# find_meta_clonotypes  (exported)
# ---------------------------------------------------------------------------

#' Identify quasi-public meta-clonotypes
#'
#' Finds TCR sequences that are shared across multiple subjects (individuals)
#' based on TCRdist neighborhoods. A meta-clonotype is defined by a center
#' TCR and a radius: all TCRs within the radius are considered part of the
#' meta-clonotype.
#'
#' @param tcr_df Data.frame with TCR columns (\code{va}, \code{vb},
#'   \code{cdr3a}, \code{cdr3b}) and a subject column.
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#' @param radius Numeric or integer vector. TCRdist radius for neighborhood
#'   definition. If a single value, used for all TCRs. If \code{NULL},
#'   per-TCR radii are computed from background ECDF.
#' @param background_df Data.frame. Background TCRs for radius computation.
#'   Required if \code{radius} is \code{NULL}.
#' @param ctrl_bkgd Numeric. Background proportion threshold for automatic
#'   radius. Default \code{1e-5}.
#' @param max_radius Numeric. Maximum radius to consider. Default \code{50L}.
#' @param min_nsubject Integer. Minimum number of distinct subjects in a
#'   neighborhood. Default \code{2L}.
#' @param subject_col Character string. Column name for subject IDs.
#'   Default \code{"subject"}.
#'
#' @return A data.frame with one row per meta-clonotype center:
#'   \describe{
#'     \item{\code{center_index}}{Row index in \code{tcr_df}.}
#'     \item{\code{va, cdr3a, vb, cdr3b}}{Center TCR sequence.}
#'     \item{\code{radius}}{TCRdist radius used.}
#'     \item{\code{K_neighbors}}{Total neighbors within radius.}
#'     \item{\code{nsubject}}{Number of distinct subjects.}
#'     \item{\code{neighbor_indices}}{Comma-separated neighbor row indices.}
#'   }
#'
#' @examples
#' \dontrun{
#' meta <- find_meta_clonotypes(tcr_df, "human", radius = 30,
#'                               subject_col = "subject")
#' }
#'
#' @seealso \code{\link{summarize_meta_clonotype}}, \code{\link{find_clumping}}
#' @export
find_meta_clonotypes <- function(tcr_df, organism,
                                  radius = NULL,
                                  background_df = NULL,
                                  ctrl_bkgd = 1e-5,
                                  max_radius = 50L,
                                  min_nsubject = 2L,
                                  subject_col = "subject") {
    n <- nrow(tcr_df)
    stopifnot(n >= 1L)

    if (!subject_col %in% colnames(tcr_df)) {
        stop(sprintf(
            "find_meta_clonotypes: column '%s' not found in tcr_df",
            subject_col), call. = FALSE)
    }

    subjects <- tcr_df[[subject_col]]

    # Compute per-TCR radii if not provided
    if (is.null(radius)) {
        if (is.null(background_df)) {
            stop("find_meta_clonotypes: 'background_df' is required when ",
                 "'radius' is NULL", call. = FALSE)
        }
        radii <- .calc_radii(tcr_df, organism, background_df,
                              ctrl_bkgd, max_radius)
    } else if (length(radius) == 1L) {
        radii <- rep(as.integer(radius), n)
    } else {
        stopifnot(length(radius) == n)
        radii <- as.integer(radius)
    }

    # Compute dense distance matrix
    dist_mat <- tcrdist_matrix(tcr_df, organism)

    # Find neighborhoods and count subjects
    results <- vector("list", n)
    neighbor_sets <- vector("list", n)

    for (i in seq_len(n)) {
        nbrs <- which(dist_mat[i, ] <= radii[i])
        nbr_subjects <- unique(subjects[nbrs])
        nsubject <- length(nbr_subjects)

        neighbor_sets[[i]] <- sort(nbrs)

        if (nsubject >= min_nsubject) {
            results[[i]] <- data.frame(
                center_index = i,
                va = tcr_df$va[i],
                cdr3a = tcr_df$cdr3a[i],
                vb = tcr_df$vb[i],
                cdr3b = tcr_df$cdr3b[i],
                radius = radii[i],
                K_neighbors = length(nbrs),
                nsubject = nsubject,
                neighbor_indices = paste(nbrs, collapse = ","),
                stringsAsFactors = FALSE
            )
        }
    }

    # Filter non-null results
    results <- results[!vapply(results, is.null, logical(1L))]

    if (length(results) == 0L) {
        return(data.frame(
            center_index = integer(0), va = character(0),
            cdr3a = character(0), vb = character(0),
            cdr3b = character(0), radius = integer(0),
            K_neighbors = integer(0), nsubject = integer(0),
            neighbor_indices = character(0),
            stringsAsFactors = FALSE
        ))
    }

    result_df <- do.call(rbind, results)
    rownames(result_df) <- NULL

    # De-duplicate: remove centers whose neighbor sets are strict subsets
    # of another center's neighbor set
    keep <- rep(TRUE, nrow(result_df))
    center_indices <- result_df$center_index

    for (a in seq_len(nrow(result_df))) {
        if (!keep[a]) next
        set_a <- neighbor_sets[[center_indices[a]]]
        for (b in seq_len(nrow(result_df))) {
            if (a == b || !keep[b]) next
            set_b <- neighbor_sets[[center_indices[b]]]
            # If a's set is a strict subset of b's set, remove a
            if (length(set_a) < length(set_b) &&
                all(set_a %in% set_b)) {
                keep[a] <- FALSE
                break
            }
        }
    }

    result_df <- result_df[keep, , drop = FALSE]
    rownames(result_df) <- NULL
    result_df
}


# ---------------------------------------------------------------------------
# summarize_meta_clonotype  (exported)
# ---------------------------------------------------------------------------

#' Summarize a meta-clonotype's composition
#'
#' Reports the V/J gene usage, CDR3 length distribution, and subject
#' representation within a meta-clonotype neighborhood.
#'
#' @param tcr_df Data.frame with TCR columns.
#' @param center_idx Integer. Row index of the center TCR.
#' @param neighbor_indices Integer vector. Row indices of neighbors.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{center}}{Named list of center TCR fields.}
#'     \item{\code{n_members}}{Number of members.}
#'     \item{\code{va_usage}}{Named integer vector of V-alpha gene counts.}
#'     \item{\code{vb_usage}}{Named integer vector of V-beta gene counts.}
#'     \item{\code{cdr3a_lengths}}{Integer vector of CDR3-alpha lengths.}
#'     \item{\code{cdr3b_lengths}}{Integer vector of CDR3-beta lengths.}
#'   }
#'
#' @examples
#' \dontrun{
#' meta <- find_meta_clonotypes(tcr_df, "mouse", radius = 30,
#'                               subject_col = "subject")
#' if (nrow(meta) > 0) {
#'     indices <- as.integer(strsplit(meta$neighbor_indices[1], ",")[[1]])
#'     summary <- summarize_meta_clonotype(
#'         tcr_df, meta$center_index[1], indices)
#'     summary$n_members
#'     summary$va_usage
#' }
#' }
#'
#' @seealso \code{\link{find_meta_clonotypes}}
#' @export
summarize_meta_clonotype <- function(tcr_df, center_idx, neighbor_indices) {
    stopifnot(center_idx >= 1L, center_idx <= nrow(tcr_df))
    neighbor_indices <- as.integer(neighbor_indices)

    members <- tcr_df[neighbor_indices, , drop = FALSE]

    list(
        center = list(
            va = tcr_df$va[center_idx],
            cdr3a = tcr_df$cdr3a[center_idx],
            vb = tcr_df$vb[center_idx],
            cdr3b = tcr_df$cdr3b[center_idx]
        ),
        n_members = length(neighbor_indices),
        va_usage = sort(table(members$va), decreasing = TRUE),
        vb_usage = sort(table(members$vb), decreasing = TRUE),
        cdr3a_lengths = nchar(members$cdr3a),
        cdr3b_lengths = nchar(members$cdr3b)
    )
}
