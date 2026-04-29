# Repertoire diversity metrics
#
# Ported from Python tcrdist3 diversity.py.
# Generalized Simpson's entropy, fuzzy diversity, clonality.


# ---------------------------------------------------------------------------
# tcr_diversity  (exported)
# ---------------------------------------------------------------------------

#' Generalized Simpson's diversity (Hill numbers framework)
#'
#' Computes the order-\code{r} generalized Simpson's entropy for a vector of
#' clonotype counts. For \code{order=2}, this is the classical Simpson's
#' diversity index. Higher orders give more weight to dominant clonotypes.
#'
#' The effective number of species (Hill number) is also returned:
#' \code{D = 1 / (1 - Z^(1/r))}.
#'
#' @param counts Integer vector. Clonotype counts (positive integers).
#' @param order Integer. Order of the diversity index. Default \code{2L}
#'   (standard Simpson's).
#' @param ci Logical. If \code{TRUE} (default), compute confidence interval.
#' @param alpha Numeric. Significance level for CI. Default \code{0.05}.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{entropy}}{Numeric. The diversity index Z_r, between 0 and 1.}
#'     \item{\code{effective_number}}{Numeric. Hill number (effective species).}
#'     \item{\code{ci_lower}}{Numeric. Lower CI bound (if \code{ci=TRUE}).}
#'     \item{\code{ci_upper}}{Numeric. Upper CI bound (if \code{ci=TRUE}).}
#'     \item{\code{order}}{Integer. The order used.}
#'   }
#'
#' @examples
#' # Uniform distribution: maximum diversity
#' tcr_diversity(rep(10, 5))
#'
#' # Single dominant clonotype: low diversity
#' tcr_diversity(c(100, 1, 1, 1))
#'
#' @seealso \code{\link{tcr_fuzzy_diversity}}, \code{\link{tcr_richness}}, \code{\link{tcr_clonality}}
#' @export
tcr_diversity <- function(counts, order = 2L, ci = TRUE, alpha = 0.05) {
    counts <- as.integer(counts)
    counts <- counts[counts > 0L]
    stopifnot(length(counts) >= 1L, order >= 2L)

    n <- sum(counts)
    S <- length(counts)
    p <- counts / n
    r <- as.integer(order)

    if (S == 1L || n < r) {
        return(list(
            entropy = 0,
            effective_number = 1,
            ci_lower = if (ci) 0 else NULL,
            ci_upper = if (ci) 0 else NULL,
            order = r
        ))
    }

    # Z_r = sum_i p_i * prod_{k=1}^{r-1} (n - c_i) / (n - k) ... no
    # Actually: Z_r = sum_i p_i * prod_{k=1}^{r-1} (c_i - 1) * ... / (n - k)
    # Corrected formula from tcrdist3:
    # Z_r = sum_i (c_i / n) * prod_{k=1}^{r-1} (1 - (c_i - 1)/(n - k))
    # Wait, let me re-derive. From the source:
    # I(p) = (1-p)^r  interpretation
    # Z_r = Σ p_i * Π_{k=1}^{r-1} ( (c_i - 1) / (n - k) )
    # This is actually the probability that r random draws are all the same species
    # Simpson's index (r=2): Z = Σ p_i * (c_i - 1)/(n - 1)
    # For higher orders: Z = Σ (c_i/n) * Π_{k=1}^{r-1} (c_i - k)/(n - k)

    # Compute Z_r
    z <- 0
    for (i in seq_along(counts)) {
        ci_val <- counts[i]
        term <- ci_val / n
        for (k in seq_len(r - 1L)) {
            if (n - k == 0) {
                term <- 0
                break
            }
            term <- term * (ci_val - k) / (n - k)
        }
        z <- z + max(0, term)
    }

    # Diversity = 1 - Z (probability of drawing r different species)
    diversity <- 1 - z

    # Effective number (Hill number)
    eff_num <- if (z > 0 && z < 1) {
        1 / z^(1 / (r - 1))
    } else if (z == 0) {
        Inf
    } else {
        1
    }

    result <- list(
        entropy = diversity,
        effective_number = eff_num,
        order = r
    )

    if (ci) {
        # Delta method CI approximation
        # Variance of Z via multinomial covariance
        # Var(Z) ≈ sum_i (dZ/dp_i)^2 * p_i*(1-p_i)/n
        # For order 2: Z = sum p_i^2, so dZ/dp_i = 2*p_i
        # Simplified: se ≈ sqrt(sum(4 * p^2 * p * (1-p)) / n)
        if (r == 2L) {
            var_z <- sum(4 * p^2 * p * (1 - p)) / n
        } else {
            # General case: numerical gradient
            var_z <- 0
            for (i in seq_along(p)) {
                term_i <- p[i]^(r - 1)
                for (k in seq_len(r - 1L)) {
                    if (n - k > 0) {
                        term_i <- term_i * (counts[i] - k) / (n - k)
                    }
                }
                grad_i <- r * term_i
                var_z <- var_z + grad_i^2 * p[i] * (1 - p[i]) / n
            }
        }
        se <- sqrt(max(0, var_z))
        z_crit <- stats::qnorm(1 - alpha / 2)
        result$ci_lower <- max(0, diversity - z_crit * se)
        result$ci_upper <- min(1, diversity + z_crit * se)
    }

    result
}


# ---------------------------------------------------------------------------
# tcr_fuzzy_diversity  (exported)
# ---------------------------------------------------------------------------

#' TCR-aware fuzzy diversity
#'
#' Computes diversity accounting for sequence similarity: two clonotypes
#' are considered "the same" if their TCRdist is within \code{threshold}.
#' This gives lower diversity for repertoires with many similar sequences.
#'
#' @param tcr_df Data.frame with TCR columns (\code{va}, \code{vb},
#'   \code{cdr3a}, \code{cdr3b}).
#' @param organism Character string (\code{"human"} or \code{"mouse"}).
#' @param threshold Numeric. Distance threshold for considering two TCRs
#'   as similar. Default \code{50}.
#' @param order Integer. Diversity order. Default \code{2L}.
#' @param counts Integer vector. Clonotype counts (one per row of
#'   \code{tcr_df}). If \code{NULL} (default), all counts are 1.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{fuzzy_diversity}}{Numeric. Fuzzy diversity, between 0 and 1.}
#'     \item{\code{standard_diversity}}{Numeric. Standard Simpson's diversity
#'       for comparison.}
#'   }
#'
#' @examples
#' \dontrun{
#' tcr_fuzzy_diversity(tcr_df, "human", threshold = 50)
#' }
#'
#' @seealso \code{\link{tcr_diversity}}, \code{\link{tcrdist_matrix}}
#' @export
tcr_fuzzy_diversity <- function(tcr_df, organism, threshold = 50,
                                 order = 2L, counts = NULL) {
    n <- nrow(tcr_df)
    if (is.null(counts)) {
        counts <- rep(1L, n)
    }
    counts <- as.integer(counts)
    stopifnot(length(counts) == n, all(counts > 0L))

    total <- sum(counts)

    # Compute sparse distance matrix within threshold
    dist_sparse <- tcrdist_sparse(tcr_df, organism, threshold = threshold)

    if (order == 2L) {
        # Analytical: fuzzy_Z = sum_{i,j} c_i * c_j * I(d(i,j) <= threshold) / total^2
        fuzzy_z <- 0
        # Diagonal contributions (always within threshold)
        fuzzy_z <- fuzzy_z + sum(as.numeric(counts)^2)

        # Off-diagonal contributions from sparse matrix
        if (inherits(dist_sparse, "dgCMatrix") || inherits(dist_sparse, "dgTMatrix")) {
            sm <- methods::as(dist_sparse, "TsparseMatrix")
            # Each entry (i,j) with i != j contributes c_i * c_j * 2 (symmetric)
            for (idx in seq_along(sm@i)) {
                i <- sm@i[idx] + 1L  # 0-based to 1-based
                j <- sm@j[idx] + 1L
                if (i < j) {  # upper triangle only, multiply by 2
                    fuzzy_z <- fuzzy_z + 2 * counts[i] * counts[j]
                }
            }
        }

        fuzzy_z <- fuzzy_z / total^2
        fuzzy_diversity <- 1 - fuzzy_z
    } else {
        # Sampling-based for higher orders
        # Sample r random indices with probability proportional to counts
        n_samples <- 10000L
        fuzzy_z_sum <- 0L

        prob <- counts / total
        for (s in seq_len(n_samples)) {
            idx <- sample(n, order, replace = TRUE, prob = prob)
            all_close <- TRUE
            for (a in seq_len(order - 1L)) {
                for (b in (a + 1L):order) {
                    i <- idx[a]
                    j <- idx[b]
                    if (i == j) next
                    d <- dist_sparse[i, j]
                    if (is.na(d) || d == 0 || d > threshold) {
                        # d==0 from sparse means > threshold (stored as 0 = not present)
                        # Need to check: sparse stores distances <= threshold
                        # Non-stored entries are > threshold
                        # Actually sparse matrix stores the distance value if <= threshold
                        # If not stored, it's 0 in the matrix but means > threshold
                        # We need to distinguish stored-0 (identical) from not-stored (>threshold)
                        # For dgCMatrix, 0 means either not stored or stored as 0
                        # Since i != j and identical TCRs have d=0 which IS stored...
                        # Actually the sparse format doesn't store explicit zeros.
                        # So d==0 for i!=j means > threshold.
                        if (i != j && d == 0) {
                            all_close <- FALSE
                            break
                        }
                    }
                }
                if (!all_close) break
            }
            if (all_close) fuzzy_z_sum <- fuzzy_z_sum + 1L
        }
        fuzzy_z <- fuzzy_z_sum / n_samples
        fuzzy_diversity <- 1 - fuzzy_z
    }

    standard <- tcr_diversity(counts, order = order, ci = FALSE)

    list(
        fuzzy_diversity = fuzzy_diversity,
        standard_diversity = standard$entropy
    )
}


# ---------------------------------------------------------------------------
# tcr_richness  (exported)
# ---------------------------------------------------------------------------

#' TCR repertoire richness
#'
#' Returns the number of unique clonotypes (species richness).
#'
#' @param counts Integer vector. Clonotype counts (positive integers).
#' @return An integer: the number of unique clonotypes.
#'
#' @examples
#' tcr_richness(c(10, 5, 3, 1, 1))  # 5
#'
#' @seealso \code{\link{tcr_diversity}}, \code{\link{tcr_clonality}}
#' @export
tcr_richness <- function(counts) {
    counts <- as.integer(counts)
    sum(counts > 0L)
}


# ---------------------------------------------------------------------------
# tcr_clonality  (exported)
# ---------------------------------------------------------------------------

#' TCR repertoire clonality
#'
#' Computes clonality as 1 minus the normalized Shannon entropy:
#' \code{1 - H / log(S)} where \code{H = -sum(p * log(p))} and \code{S}
#' is the number of unique clonotypes. Values near 0 indicate a uniform
#' (diverse) repertoire; values near 1 indicate a dominated repertoire.
#'
#' @param counts Integer vector. Clonotype counts (positive integers).
#' @return A numeric scalar between 0 and 1.
#'
#' @examples
#' # Uniform: clonality near 0
#' tcr_clonality(rep(10, 5))
#'
#' # Dominated: clonality near 1
#' tcr_clonality(c(1000, 1, 1))
#'
#' @seealso \code{\link{tcr_diversity}}, \code{\link{tcr_richness}}
#' @export
tcr_clonality <- function(counts) {
    counts <- as.integer(counts)
    counts <- counts[counts > 0L]
    S <- length(counts)
    if (S <= 1L) return(1)

    p <- counts / sum(counts)
    H <- -sum(p * log(p))
    1 - H / log(S)
}
