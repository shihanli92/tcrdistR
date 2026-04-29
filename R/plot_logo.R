# CDR3 sequence logos, junction bars, and gene usage plots
#
# Ported from rconga R/plot_logo.R (Python CoNGA make_tcr_logo.py).
# Dependencies: ggplot2 (Suggests), ggseqlogo (Suggests), patchwork (Suggests)


# ---------------------------------------------------------------------------
# .align_cdr3_regions  (internal)
# ---------------------------------------------------------------------------

#' Align two CDR3 amino acid sequences by inserting gaps in the shorter one
#'
#' Inserts gap characters into the shorter sequence at the position that
#' maximises the BLOSUM62 alignment score against the longer sequence. When
#' both sequences are the same length, returns them unchanged.
#'
#' @param a Character string. First CDR3 amino acid sequence.
#' @param b Character string. Second CDR3 amino acid sequence.
#' @param gap_character Character string of length 1. Gap character to insert.
#'   Default \code{"-"}.
#' @return A named list with elements \code{a} and \code{b}, each a character
#'   string of equal length (the aligned sequences).
#' @keywords internal
#' @noRd
.align_cdr3_regions <- function(a, b, gap_character = "-") {
    a_chars <- strsplit(a, "", fixed = TRUE)[[1L]]
    b_chars <- strsplit(b, "", fixed = TRUE)[[1L]]
    len_a <- length(a_chars)
    len_b <- length(b_chars)

    if (len_a == len_b) {
        return(list(a = a, b = b))
    }

    if (len_a < len_b) {
        s0 <- a_chars
        s1 <- b_chars
    } else {
        s0 <- b_chars
        s1 <- a_chars
    }

    len_s0 <- length(s0)
    len_s1 <- length(s1)
    lendiff <- len_s1 - len_s0

    best_score <- -1000
    best_gappos <- 0L

    valid_aa <- rownames(BLOSUM62)

    if (len_s0 >= 2L) {
        for (gappos in seq.int(0L, len_s0 - 2L)) {
            score <- 0
            for (i in seq.int(0L, gappos)) {
                aa_s0 <- s0[i + 1L]
                aa_s1 <- s1[i + 1L]
                if (aa_s0 %in% valid_aa && aa_s1 %in% valid_aa) {
                    score <- score + BLOSUM62[aa_s0, aa_s1]
                }
            }
            if (gappos + 1L <= len_s0 - 1L) {
                for (i in seq.int(gappos + 1L, len_s0 - 1L)) {
                    aa_s0 <- s0[i + 1L]
                    aa_s1 <- s1[i + lendiff + 1L]
                    if (aa_s0 %in% valid_aa && aa_s1 %in% valid_aa) {
                        score <- score + BLOSUM62[aa_s0, aa_s1]
                    }
                }
            }
            if (score > best_score) {
                best_score <- score
                best_gappos <- gappos
            }
        }
    }

    left_part <- paste(s0[seq_len(best_gappos + 1L)], collapse = "")
    gap_part <- paste(rep(gap_character, lendiff), collapse = "")
    right_part <- if (best_gappos + 2L <= len_s0) {
        paste(s0[seq.int(best_gappos + 2L, len_s0)], collapse = "")
    } else {
        ""
    }
    s0_aligned <- paste0(left_part, gap_part, right_part)

    if (len_a < len_b) {
        list(a = s0_aligned, b = paste(s1, collapse = ""))
    } else {
        list(a = paste(s1, collapse = ""), b = s0_aligned)
    }
}


# ---------------------------------------------------------------------------
# .build_cdr3_pwm  (internal)
# ---------------------------------------------------------------------------

#' Build a Position Weight Matrix from CDR3 amino acid sequences
#'
#' Aligns a set of CDR3 sequences to a center sequence and computes
#' per-position amino acid frequencies. The center is chosen as the sequence
#' with the highest total pairwise BLOSUM62 score among sequences of the most
#' common length.
#'
#' @param cdr3_seqs Character vector. CDR3 amino acid sequences.
#' @param trim Logical. If \code{TRUE}, remove the first 3 and last 2
#'   characters from each CDR3 (conserved C...F/W residues) before alignment.
#' @param gap_character Character string of length 1. Gap character.
#' @param nucseq_src List of character vectors for junction PWM (optional).
#' @return A named list with \code{pwm}, \code{aligned_seqs},
#'   \code{center_idx}, and optionally \code{junction_pwm}.
#' @keywords internal
#' @noRd
.build_cdr3_pwm <- function(cdr3_seqs, trim = TRUE, gap_character = "-",
                             nucseq_src = NULL) {
    stopifnot(is.character(cdr3_seqs), length(cdr3_seqs) >= 1L)

    build_junction <- !is.null(nucseq_src)

    # ---- Trim conserved residues ----
    if (trim) {
        cdr3_seqs <- vapply(cdr3_seqs, function(s) {
            nc <- nchar(s)
            if (nc <= 5L) return(s)
            substr(s, 4L, nc - 2L)
        }, character(1L), USE.NAMES = FALSE)

        if (build_junction) {
            nucseq_src <- lapply(nucseq_src, function(src) {
                if (is.null(src)) return(NULL)
                n_nuc <- length(src)
                if (n_nuc <= 15L) return(src)
                src[10L:(n_nuc - 6L)]
            })
        }
    }

    n_seqs <- length(cdr3_seqs)

    # ---- Single sequence: trivial PWM ----
    if (n_seqs == 1L) {
        chars <- strsplit(cdr3_seqs, "", fixed = TRUE)[[1L]]
        L <- length(chars)
        aa_labels <- c(AMINO_ACIDS, gap_character)
        pwm <- matrix(0, nrow = length(aa_labels), ncol = L,
                       dimnames = list(aa_labels, NULL))
        for (j in seq_len(L)) {
            if (chars[j] %in% aa_labels) {
                pwm[chars[j], j] <- 1.0
            }
        }
        result <- list(pwm = pwm, aligned_seqs = cdr3_seqs, center_idx = 1L)
        if (build_junction && !is.null(nucseq_src[[1L]])) {
            src <- nucseq_src[[1L]]
            junc_labels <- unique(c("V", "N", "N1", "D", "N2", "J",
                                    gap_character))
            jpwm <- matrix(0, nrow = length(junc_labels), ncol = 3L * L,
                           dimnames = list(junc_labels, NULL))
            for (k in seq_along(src)) {
                if (src[k] %in% junc_labels) {
                    jpwm[src[k], k] <- 1.0
                }
            }
            result$junction_pwm <- jpwm
        }
        return(result)
    }

    # ---- Find center sequence ----
    seq_lengths <- nchar(cdr3_seqs)
    length_table <- table(seq_lengths)
    most_common_length <- as.integer(
        names(length_table)[which.max(length_table)])
    candidate_idx <- which(seq_lengths == most_common_length)

    if (length(candidate_idx) == 1L) {
        center_idx <- candidate_idx[1L]
    } else {
        best_total <- -Inf
        center_idx <- candidate_idx[1L]

        for (ci in candidate_idx) {
            total_score <- 0
            for (oi in seq_len(n_seqs)) {
                if (oi == ci) next
                aligned <- .align_cdr3_regions(cdr3_seqs[ci], cdr3_seqs[oi],
                                               gap_character)
                a_chars <- strsplit(aligned$a, "", fixed = TRUE)[[1L]]
                b_chars <- strsplit(aligned$b, "", fixed = TRUE)[[1L]]
                for (j in seq_along(a_chars)) {
                    aa <- a_chars[j]
                    bb <- b_chars[j]
                    if (aa != gap_character && bb != gap_character &&
                        aa %in% AMINO_ACIDS && bb %in% AMINO_ACIDS) {
                        total_score <- total_score + BLOSUM62[aa, bb]
                    }
                }
            }
            if (total_score > best_total) {
                best_total <- total_score
                center_idx <- ci
            }
        }
    }

    center_cdr3 <- cdr3_seqs[center_idx]
    center_len <- nchar(center_cdr3)

    # ---- Align all sequences to center and build PWM ----
    aa_labels <- c(AMINO_ACIDS, gap_character)
    pwm <- matrix(0, nrow = length(aa_labels), ncol = center_len,
                   dimnames = list(aa_labels, NULL))

    junction_pwm <- NULL
    if (build_junction) {
        junc_labels <- c("V", "N", "N1", "D", "N2", "J", gap_character)
        junction_pwm <- matrix(0, nrow = length(junc_labels),
                               ncol = 3L * center_len,
                               dimnames = list(junc_labels, NULL))
    }

    aligned_seqs <- character(n_seqs)

    for (si in seq_len(n_seqs)) {
        aligned <- .align_cdr3_regions(center_cdr3, cdr3_seqs[si],
                                       gap_character)
        a_chars <- strsplit(aligned$a, "", fixed = TRUE)[[1L]]
        b_chars <- strsplit(aligned$b, "", fixed = TRUE)[[1L]]

        member_junc <- if (build_junction) nucseq_src[[si]] else NULL
        b_is_gap <- (b_chars == gap_character)
        b_cum_gaps <- cumsum(b_is_gap)

        gaps_in_a <- 0L
        kept_chars <- character(center_len)
        for (j in seq_along(a_chars)) {
            if (a_chars[j] == gap_character) {
                gaps_in_a <- gaps_in_a + 1L
                next
            }
            pwmpos <- j - gaps_in_a
            b_aa <- b_chars[j]
            kept_chars[pwmpos] <- b_aa
            if (b_aa %in% aa_labels) {
                pwm[b_aa, pwmpos] <- pwm[b_aa, pwmpos] + 1
            }

            if (!is.null(member_junc)) {
                for (k in 0L:2L) {
                    jpwm_col <- 3L * (pwmpos - 1L) + k + 1L
                    if (b_aa == gap_character) {
                        junction_pwm[gap_character, jpwm_col] <-
                            junction_pwm[gap_character, jpwm_col] + 1
                    } else {
                        bpos_0 <- (j - 1L) - b_cum_gaps[j]
                        nuc_idx <- 3L * bpos_0 + k + 1L
                        if (nuc_idx >= 1L &&
                            nuc_idx <= length(member_junc)) {
                            bsrc <- member_junc[nuc_idx]
                            if (bsrc %in% rownames(junction_pwm)) {
                                junction_pwm[bsrc, jpwm_col] <-
                                    junction_pwm[bsrc, jpwm_col] + 1
                            }
                        }
                    }
                }
            }
        }
        aligned_seqs[si] <- paste(kept_chars, collapse = "")
    }

    # ---- Normalize columns ----
    col_sums <- colSums(pwm)
    for (j in seq_len(ncol(pwm))) {
        if (col_sums[j] > 0) {
            pwm[, j] <- pwm[, j] / col_sums[j]
        }
    }

    result <- list(pwm = pwm, aligned_seqs = aligned_seqs,
                   center_idx = center_idx)

    if (build_junction && !is.null(junction_pwm)) {
        junc_col_sums <- colSums(junction_pwm)
        for (j in seq_len(ncol(junction_pwm))) {
            if (junc_col_sums[j] > 0) {
                junction_pwm[, j] <- junction_pwm[, j] / junc_col_sums[j]
            }
        }
        result$junction_pwm <- junction_pwm
    }

    result
}


# ---------------------------------------------------------------------------
# plot_cdr3_logo  (exported)
# ---------------------------------------------------------------------------

#' Plot a CDR3 sequence logo
#'
#' Renders a sequence logo from CDR3 amino acid sequences using
#' \pkg{ggseqlogo}. Two methods are available: \code{"bits"} shows information
#' content per position, \code{"prob"} shows amino acid frequencies.
#'
#' @param cdr3_seqs Character vector. CDR3 amino acid sequences.
#' @param chain Character string. \code{"alpha"} or \code{"beta"}.
#' @param trim Logical. If \code{TRUE}, trim conserved C/F residues before
#'   alignment (first 3, last 2 characters). Default \code{TRUE}.
#' @param method Character string. \code{"prob"} for frequency-based logo,
#'   \code{"bits"} for information content. Default \code{"prob"}.
#' @param gap_character Character string. Gap character for alignment.
#'   Default \code{"-"}.
#' @param title Optional plot title. If \code{NULL}, defaults to the chain
#'   name with Greek letter.
#' @param show_junction_bars Logical. If \code{TRUE} and \code{nucseq_src}
#'   is provided, stack junction bars below the logo via \pkg{patchwork}.
#' @param nucseq_src List of character vectors with nucleotide source labels
#'   (one per CDR3). Required for junction bars.
#'
#' @return A \code{ggplot} object (or \code{patchwork} object if junction bars
#'   are included).
#'
#' @examples
#' \dontrun{
#' seqs <- c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF",
#'           "CAVRDSSYKLIF", "CAVKDSYKLIF")
#' plot_cdr3_logo(seqs, chain = "alpha", method = "bits")
#' plot_cdr3_logo(seqs, chain = "beta", method = "prob")
#' }
#'
#' @seealso \code{\link{plot_junction_bars}}, \code{\link{plot_gene_usage}}
#' @export
plot_cdr3_logo <- function(cdr3_seqs,
                            chain = c("alpha", "beta"),
                            trim = TRUE,
                            method = c("prob", "bits"),
                            gap_character = "-",
                            title = NULL,
                            nucseq_src = NULL,
                            show_junction_bars = FALSE) {
    .check_ggplot2("plot_cdr3_logo()")
    if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
        stop("Package 'ggseqlogo' is required for plot_cdr3_logo(). ",
             "Install it with: install.packages(\"ggseqlogo\")",
             call. = FALSE)
    }

    chain <- match.arg(chain)
    method <- match.arg(method)
    stopifnot(is.character(cdr3_seqs), length(cdr3_seqs) >= 1L)

    pwm_result <- .build_cdr3_pwm(cdr3_seqs, trim = trim,
                                   gap_character = gap_character,
                                   nucseq_src = nucseq_src)

    if (is.null(title)) {
        chain_label <- if (chain == "alpha") "CDR3\u03b1" else "CDR3\u03b2"
        title <- chain_label
    }

    if (method == "bits") {
        aligned <- pwm_result$aligned_seqs
        if (gap_character != "-") {
            aligned <- gsub(gap_character, "-", aligned, fixed = TRUE)
        }
        p <- ggseqlogo::ggseqlogo(aligned, method = "bits",
                                   seq_type = "aa") +
            ggplot2::ggtitle(title)
    } else {
        pwm <- pwm_result$pwm
        gap_row <- which(rownames(pwm) == gap_character)
        if (length(gap_row) > 0L) {
            pwm <- pwm[-gap_row, , drop = FALSE]
        }
        p <- ggseqlogo::ggseqlogo(pwm, method = "custom",
                                   seq_type = "aa") +
            ggplot2::ggtitle(title)
    }

    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
            axis.text.x = ggplot2::element_text(size = 7),
            axis.text.y = ggplot2::element_text(size = 7),
            legend.position = "none"
        )

    if (show_junction_bars && !is.null(pwm_result$junction_pwm)) {
        if (!requireNamespace("patchwork", quietly = TRUE)) {
            stop("Package 'patchwork' is required for junction bars.",
                 call. = FALSE)
        }
        jbar <- plot_junction_bars(pwm_result$junction_pwm,
                                    chain = chain,
                                    gap_character = gap_character)
        p <- patchwork::wrap_plots(p, jbar, ncol = 1L, heights = c(3, 1))
    }

    p
}


# ---------------------------------------------------------------------------
# plot_junction_bars  (exported)
# ---------------------------------------------------------------------------

#' Plot junction bars showing V/N/D/J nucleotide composition
#'
#' Renders stacked coloured bars showing the fraction of nucleotides derived
#' from V-gene, N-region, D-gene, or J-gene segments at each codon position
#' within a CDR3 alignment. Designed to sit below a CDR3 sequence logo.
#'
#' @param junction_pwm Numeric matrix. Rows are junction source labels,
#'   columns are nucleotide positions (3 per codon).
#' @param chain Character string. \code{"alpha"} or \code{"beta"}.
#' @param gap_character Character string. Gap character. Default \code{"-"}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{plot_cdr3_logo}}
#' @export
plot_junction_bars <- function(junction_pwm,
                                chain = c("alpha", "beta"),
                                gap_character = "-") {
    .check_ggplot2("plot_junction_bars()")

    chain <- match.arg(chain)

    if (chain == "alpha") {
        src_order <- c("V", "N", "J")
    } else {
        src_order <- c("V", "N1", "D", "N2", "J")
    }
    src_colors <- c(
        "V"  = "#C0C0C0",
        "N"  = "red",
        "N1" = "red",
        "D"  = "black",
        "N2" = "red",
        "J"  = "dimgray"
    )

    n_cols <- ncol(junction_pwm)

    df_list <- vector("list", length(src_order) * n_cols)
    idx <- 0L
    for (col_i in seq_len(n_cols)) {
        for (src in src_order) {
            idx <- idx + 1L
            frac <- if (src %in% rownames(junction_pwm)) {
                junction_pwm[src, col_i]
            } else {
                0
            }
            df_list[[idx]] <- data.frame(
                pos = col_i,
                source = src,
                fraction = frac,
                stringsAsFactors = FALSE
            )
        }
    }
    df <- do.call(rbind, df_list)
    df$source <- factor(df$source, levels = src_order)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$pos, y = .data$fraction,
                                      fill = .data$source)) +
        ggplot2::geom_col(position = "stack", width = 0.9) +
        ggplot2::scale_fill_manual(
            values = src_colors[src_order],
            breaks = src_order,
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            breaks = NULL,
            expand = ggplot2::expansion(mult = 0.01)
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = 0)
        ) +
        ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            legend.position = "none",
            plot.margin = ggplot2::margin(0, 2, 2, 2)
        )
}


# ---------------------------------------------------------------------------
# plot_gene_usage  (exported)
# ---------------------------------------------------------------------------

#' Plot V-gene or J-gene usage frequencies
#'
#' Horizontal bar chart of gene frequencies from a TCR data.frame column.
#' Optionally strips allele suffixes (e.g., \code{"TRAV1-2*01"} becomes
#' \code{"TRAV1-2"}).
#'
#' @param tcr_df Data.frame containing TCR data.
#' @param gene_col Character string. Column name to tally (e.g., \code{"va"},
#'   \code{"vb"}, \code{"ja"}).
#' @param strip_allele Logical. If \code{TRUE} (default), strip allele
#'   suffixes (everything after \code{"*"}) before tallying.
#' @param max_genes Integer. Maximum number of genes to display. Default
#'   \code{20L}.
#' @param title Optional plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' tcr_df <- data.frame(va = c("TRAV1-2*01", "TRAV1-2*02", "TRAV12-1*01",
#'                              "TRAV1-2*01", "TRAV12-1*01"))
#' plot_gene_usage(tcr_df, "va", title = "V-alpha usage")
#' }
#'
#' @seealso \code{\link{plot_cdr3_logo}}, \code{\link{plot_tcrdist_dendrogram}}
#' @export
plot_gene_usage <- function(tcr_df, gene_col, strip_allele = TRUE,
                             max_genes = 20L, title = NULL) {
    .check_ggplot2("plot_gene_usage()")

    stopifnot(gene_col %in% colnames(tcr_df))

    genes <- tcr_df[[gene_col]]
    genes <- genes[!is.na(genes) & nzchar(genes)]

    if (length(genes) == 0L) {
        return(ggplot2::ggplot() +
                   ggplot2::theme_void() +
                   ggplot2::ggtitle(title %||% "No gene data"))
    }

    if (strip_allele) {
        genes <- sub("\\*.*$", "", genes)
    }

    freq_table <- sort(table(genes), decreasing = TRUE)
    if (length(freq_table) > max_genes) {
        freq_table <- freq_table[seq_len(max_genes)]
    }

    df <- data.frame(
        gene = factor(names(freq_table), levels = rev(names(freq_table))),
        count = as.integer(freq_table),
        stringsAsFactors = FALSE
    )

    if (is.null(title)) {
        title <- paste0(gene_col, " usage")
    }

    ggplot2::ggplot(df, ggplot2::aes(x = .data$gene, y = .data$count)) +
        ggplot2::geom_col(fill = "#1f77b4", width = 0.7) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Count", title = title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
        )
}
