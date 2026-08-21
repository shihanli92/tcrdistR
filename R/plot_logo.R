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
#' @param return_junction_pwm Logical. If \code{TRUE}, return a list with
#'   components \code{plot} (the ggplot logo) and \code{junction_pwm} (the
#'   numeric matrix, or \code{NULL}). Default \code{FALSE}.
#'
#' @return A \code{ggplot} object (or \code{patchwork} object if junction bars
#'   are included). When \code{return_junction_pwm = TRUE}, a list with
#'   \code{plot} and \code{junction_pwm}.
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
                            show_junction_bars = FALSE,
                            return_junction_pwm = FALSE) {
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

    if (return_junction_pwm) {
        # Override x limits to align with junction bars
        n_aa <- ncol(pwm_result$pwm)
        p <- p + ggplot2::scale_x_continuous(
            limits = c(0.5, n_aa + 0.5), expand = c(0, 0))
        return(list(plot = p, junction_pwm = pwm_result$junction_pwm))
    }

    if (show_junction_bars && !is.null(pwm_result$junction_pwm)) {
        if (!requireNamespace("patchwork", quietly = TRUE)) {
            stop("Package 'patchwork' is required for junction bars.",
                 call. = FALSE)
        }
        # Override logo x limits to match junction bar scale
        n_aa <- ncol(pwm_result$pwm)
        p <- p + ggplot2::scale_x_continuous(
            limits = c(0.5, n_aa + 0.5), expand = c(0, 0))
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
    n_aa   <- n_cols %/% 3L

    # Map nucleotide positions to amino-acid x-scale so bars align with
    # ggseqlogo letters (which sit at integer positions 1..n_aa).
    # Each codon's 3 nucleotides tile within the amino acid's [k-0.5, k+0.5].
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
                pos = (col_i - 0.5) / 3 + 0.5,
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
        ggplot2::geom_col(position = "stack", width = 1 / 3) +
        ggplot2::scale_fill_manual(
            values = src_colors[src_order],
            breaks = src_order,
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            breaks = NULL,
            limits = c(0.5, n_aa + 0.5),
            expand = c(0, 0)
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
#' @param tcr_rep A \code{\linkS4class{TCRrep}} object.  If provided, the
#'   clone data.frame is used as \code{tcr_df}.  Cannot be combined with
#'   \code{tcr_df}.
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
plot_gene_usage <- function(tcr_df = NULL, gene_col, strip_allele = TRUE,
                             max_genes = 20L, title = NULL,
                             tcr_rep = NULL) {
    .check_ggplot2("plot_gene_usage()")

    if (!is.null(tcr_rep)) {
        if (!is.null(tcr_df))
            stop("provide 'tcr_rep' or 'tcr_df', not both", call. = FALSE)
        tcr_df <- .extract_from_tcr_rep(tcr_rep)$clone_df
    }

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


# ---------------------------------------------------------------------------
# plot_gene_alluvial  (exported)
# ---------------------------------------------------------------------------

#' Alluvial plot of TCR gene segment usage
#'
#' Draws an alluvial (Sankey) diagram showing how clonotypes flow across
#' gene segments.  The default axis order is
#' \code{ja} \eqn{\rightarrow} \code{va} \eqn{\rightarrow} \code{vb}
#' \eqn{\rightarrow} \code{jb}, but any subset or reordering can be
#' specified via \code{axes}.
#'
#' Requires the \pkg{ggalluvial} package (Suggests).
#'
#' @param tcr_df Data.frame with gene columns (e.g. \code{va}, \code{ja},
#'   \code{vb}, \code{jb}).
#' @param axes Character vector of column names to use as axes (strata),
#'   in display order.  Default \code{c("ja", "va", "vb", "jb")}.
#' @param strip_allele Logical.  If \code{TRUE} (default), strip allele
#'   suffixes before tallying.
#' @param max_genes Integer.  Per-axis cap: only the top \code{max_genes}
#'   genes are shown; the rest are collapsed to \code{"Other"}.
#'   Default \code{15L}.
#' @param title Optional plot title.
#' @param color_by Column name (character string) to color the alluvial
#'   flows by.  Must be one of the \code{axes}.  Set to \code{NULL} for
#'   a uniform fill color (see \code{fill_color}).  Default \code{NULL}
#'   (uniform fill).
#' @param fill_color Single color string used for all flows when
#'   \code{color_by} is \code{NULL}.  Default \code{"black"}.  Ignored
#'   when \code{color_by} is set.
#' @param alpha Numeric.  Flow opacity.  Default \code{0.4}.
#' @param palette Character vector of colors, or \code{NULL} for the
#'   default \code{.tcrdistR_palette}.  Only used when \code{color_by}
#'   is not \code{NULL}.
#' @param facet_by Optional column name (character string) or vector of
#'   values to facet by.  When a single string that matches a column in
#'   \code{tcr_df}, that column is used for \code{facet_wrap}.  Otherwise
#'   treated as a vector of facet labels (must have length
#'   \code{nrow(tcr_df)}).
#' @param tcr_rep A \code{\linkS4class{TCRrep}} object.  If provided, the
#'   clone data.frame is used as \code{tcr_df}.  Cannot be combined with
#'   \code{tcr_df}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' data(dash)
#' plot_gene_alluvial(dash)
#' plot_gene_alluvial(dash, color_by = "va")
#' plot_gene_alluvial(dash, fill_color = "steelblue", alpha = 0.5)
#' plot_gene_alluvial(dash, facet_by = "epitope")
#' }
#'
#' @seealso \code{\link{plot_gene_usage}}, \code{\link{plot_cdr3_length}}
#' @export
plot_gene_alluvial <- function(tcr_df = NULL,
                                axes = c("ja", "va", "vb", "jb"),
                                mode = c("sankey", "alluvial"),
                                strip_allele = TRUE,
                                max_genes = 15L,
                                title = NULL,
                                color_by = NULL,
                                fill_color = "black",
                                alpha = 0.4,
                                palette = NULL,
                                facet_by = NULL,
                                tcr_rep = NULL) {
    .check_ggplot2("plot_gene_alluvial()")
    mode <- match.arg(mode)

    if (mode == "alluvial" &&
        !requireNamespace("ggalluvial", quietly = TRUE)) {
        stop("Package 'ggalluvial' is required for mode='alluvial'. ",
             "Install it with: install.packages(\"ggalluvial\")",
             call. = FALSE)
    }

    if (!is.null(tcr_rep)) {
        if (!is.null(tcr_df))
            stop("provide 'tcr_rep' or 'tcr_df', not both", call. = FALSE)
        tcr_df <- .extract_from_tcr_rep(tcr_rep)$clone_df
    }

    stopifnot(length(axes) >= 2L)
    missing_cols <- setdiff(axes, colnames(tcr_df))
    if (length(missing_cols) > 0L) {
        stop("tcr_df is missing columns: ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    # Resolve facet_by: column name lookup or raw vector
    facet_col <- NULL
    if (!is.null(facet_by)) {
        if (is.character(facet_by) && length(facet_by) == 1L &&
            facet_by %in% colnames(tcr_df)) {
            facet_col <- ".facet"
            facet_vec <- tcr_df[[facet_by]]
        } else {
            if (length(facet_by) != nrow(tcr_df)) {
                stop("'facet_by' must be a column name or a vector of length ",
                     "nrow(tcr_df)", call. = FALSE)
            }
            facet_col <- ".facet"
            facet_vec <- facet_by
        }
    }

    # Build working data.frame
    df <- tcr_df[, axes, drop = FALSE]
    if (!is.null(facet_col)) df[[facet_col]] <- facet_vec
    df <- df[stats::complete.cases(df), , drop = FALSE]

    if (nrow(df) == 0L) {
        return(ggplot2::ggplot() + ggplot2::theme_void() +
                   ggplot2::ggtitle(title %||% "No gene data"))
    }

    if (strip_allele) {
        for (col in axes) df[[col]] <- sub("\\*.*$", "", df[[col]])
    }

    for (col in axes) {
        freq <- sort(table(df[[col]]), decreasing = TRUE)
        if (length(freq) > max_genes) {
            keep <- names(freq)[seq_len(max_genes)]
            df[[col]][!df[[col]] %in% keep] <- "Other"
        }
    }

    # Validate color_by
    use_fill_aes <- !is.null(color_by)
    if (use_fill_aes && !color_by %in% axes) {
        stop("'color_by' must be one of the axes: ",
             paste(axes, collapse = ", "), call. = FALSE)
    }

    # ---- Dispatch by mode ----------------------------------------------------
    if (mode == "alluvial") {
        p <- .gene_alluvial_alluvial(df, axes, facet_col, color_by,
                                      fill_color, alpha, palette, title)
    } else {
        p <- .gene_alluvial_sankey(df, axes, facet_col, color_by,
                                    fill_color, alpha, palette, title)
    }

    if (!is.null(facet_col)) {
        p <- p + ggplot2::facet_wrap(
            ggplot2::vars(.data[[facet_col]]),
            scales = "free_y"
        )
    }

    p
}


# -- alluvial mode (ggalluvial) ------------------------------------------------

.gene_alluvial_alluvial <- function(df, axes, facet_col, color_by,
                                     fill_color, alpha, palette, title) {
    group_cols <- axes
    if (!is.null(facet_col)) group_cols <- c(group_cols, facet_col)
    agg <- stats::aggregate(
        list(Freq = rep(1L, nrow(df))),
        by = as.list(df[, group_cols, drop = FALSE]),
        FUN = sum
    )
    for (col in axes) {
        freq <- sort(tapply(agg$Freq, agg[[col]], sum), decreasing = TRUE)
        agg[[col]] <- factor(agg[[col]], levels = names(freq))
    }

    use_fill_aes <- !is.null(color_by)

    p <- ggplot2::ggplot(
        agg,
        ggplot2::aes(axis1 = .data[[axes[1L]]],
                     axis2 = .data[[axes[2L]]],
                     y = .data$Freq)
    )
    if (length(axes) >= 3L) p <- p + ggplot2::aes(axis3 = .data[[axes[3L]]])
    if (length(axes) >= 4L) p <- p + ggplot2::aes(axis4 = .data[[axes[4L]]])

    if (use_fill_aes) {
        n_levels <- length(levels(agg[[color_by]]))
        if (is.null(palette)) palette <- .tcrdistR_palette(n_levels)
        p <- p +
            ggalluvial::geom_alluvium(
                ggplot2::aes(fill = .data[[color_by]]),
                width = 1 / 4, alpha = alpha
            ) +
            ggplot2::scale_fill_manual(values = palette)
    } else {
        p <- p +
            ggalluvial::geom_alluvium(
                width = 1 / 4, fill = fill_color, alpha = alpha
            )
    }

    p +
        ggalluvial::geom_stratum(width = 1 / 4, fill = "grey90",
                                  color = "grey40") +
        ggplot2::geom_text(
            stat = ggalluvial::StatStratum,
            ggplot2::aes(label = ggplot2::after_stat(.data$stratum)),
            size = 2.5
        ) +
        ggplot2::scale_x_discrete(limits = axes,
                                   expand = c(0.15, 0.05)) +
        ggplot2::labs(y = "Count", title = title,
                      fill = if (use_fill_aes) color_by else NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        )
}


# -- sankey mode (pure ggplot2) ------------------------------------------------

.gene_alluvial_sankey <- function(df, axes, facet_col, color_by,
                                   fill_color, alpha, palette, title) {
    n_axes   <- length(axes)
    n_curves <- 50L
    stratum_w <- 0.12

    # -- compute per-axis gene order (by overall frequency, descending) --------
    axis_levels <- lapply(axes, function(col) {
        freq <- sort(table(df[[col]]), decreasing = TRUE)
        names(freq)
    })
    names(axis_levels) <- axes

    # -- build strata layout: cumulative y-positions per facet -----------------
    facet_vals <- if (!is.null(facet_col)) unique(df[[facet_col]]) else list(NULL)

    all_strata  <- list()
    all_flows   <- list()

    for (fv in facet_vals) {
        if (!is.null(facet_col)) {
            sub <- df[df[[facet_col]] == fv, , drop = FALSE]
        } else {
            sub <- df
        }
        total_n <- nrow(sub)

        # Per-axis: compute gene counts & cumulative y positions
        axis_info <- list()
        for (ai in seq_along(axes)) {
            col <- axes[ai]
            lvls <- axis_levels[[col]]
            counts <- table(factor(sub[[col]], levels = lvls))
            yend <- cumsum(as.integer(counts))
            ystart <- c(0L, yend[-length(yend)])
            axis_info[[col]] <- data.frame(
                gene   = lvls,
                ystart = ystart,
                yend   = yend,
                count  = as.integer(counts),
                stringsAsFactors = FALSE
            )
        }

        # Build strata rectangles
        for (ai in seq_along(axes)) {
            col <- axes[ai]
            info <- axis_info[[col]]
            info <- info[info$count > 0L, , drop = FALSE]
            sr <- data.frame(
                xmin  = ai - stratum_w,
                xmax  = ai + stratum_w,
                ymin  = info$ystart,
                ymax  = info$yend,
                label = info$gene,
                axis  = col,
                stringsAsFactors = FALSE
            )
            if (!is.null(facet_col)) sr[[facet_col]] <- fv
            all_strata[[length(all_strata) + 1L]] <- sr
        }

        # Build flow polygons between each adjacent pair
        for (pi in seq_len(n_axes - 1L)) {
            left_col  <- axes[pi]
            right_col <- axes[pi + 1L]
            left_info  <- axis_info[[left_col]]
            right_info <- axis_info[[right_col]]

            # Pairwise counts
            pair_tab <- table(
                factor(sub[[left_col]],  levels = left_info$gene),
                factor(sub[[right_col]], levels = right_info$gene)
            )

            # Track consumed y per stratum
            left_used  <- stats::setNames(left_info$ystart,  left_info$gene)
            right_used <- stats::setNames(right_info$ystart, right_info$gene)

            for (lg in left_info$gene) {
                for (rg in right_info$gene) {
                    cnt <- as.integer(pair_tab[lg, rg])
                    if (cnt == 0L) next

                    ly0 <- left_used[lg]
                    ly1 <- ly0 + cnt
                    ry0 <- right_used[rg]
                    ry1 <- ry0 + cnt

                    left_used[lg]  <- ly1
                    right_used[rg] <- ry1

                    # Sigmoid curve
                    t <- seq(0, 1, length.out = n_curves)
                    s <- 3 / (1 + exp(-12 * (t - 0.5)))  # scaled sigmoid
                    s <- (s - min(s)) / (max(s) - min(s))
                    xvals <- pi + stratum_w + t * (1 - 2 * stratum_w)

                    top    <- ly1 + s * (ry1 - ly1)
                    bottom <- ly0 + s * (ry0 - ly0)

                    poly <- data.frame(
                        x = c(xvals, rev(xvals)),
                        y = c(top, rev(bottom)),
                        left_gene  = lg,
                        right_gene = rg,
                        pair       = pi,
                        stringsAsFactors = FALSE
                    )
                    if (!is.null(facet_col)) poly[[facet_col]] <- fv
                    all_flows[[length(all_flows) + 1L]] <- poly
                }
            }
        }
    }

    strata_df <- do.call(rbind, all_strata)
    flow_df   <- do.call(rbind, all_flows)

    # Assign group IDs for geom_polygon
    flow_df$.gid <- rep(seq_along(all_flows),
                        vapply(all_flows, nrow, integer(1L)))

    # -- color logic -----------------------------------------------------------
    use_fill_aes <- !is.null(color_by)
    if (use_fill_aes) {
        # color_by is always the left-side gene of the first pair
        # map the left_gene column for pair==1, otherwise the flow inherits
        # from whatever stratum it originated at on the color_by axis
        color_axis_idx <- match(color_by, axes)
        # For each flow, determine the gene at the color_by axis
        # Flows have a pair index; if color_by axis == pair's left, use left_gene
        # if color_by axis == pair's right, use right_gene
        flow_df$.color <- ifelse(
            flow_df$pair == color_axis_idx,
            flow_df$left_gene,
            ifelse(flow_df$pair == color_axis_idx - 1L,
                   flow_df$right_gene,
                   NA_character_)
        )
        # For axes not adjacent to color_by, propagate via left/right chain
        # Simpler approach: just color by left_gene of pair containing color_by
        # For a clean look, we color ALL flows by their gene at the color_by axis
        # This requires tracing — for simplicity, color each pair's left gene
        # when color_by is the first axis, or right gene when it's the last
        if (color_axis_idx == 1L) {
            flow_df$.color <- flow_df$left_gene
            # Only accurate for pair 1; other pairs need propagation
            # For sankey, flows don't chain, so just color per-segment
        }
        # Simplest correct approach: color by left_gene for all flows where
        # pair >= color_axis_idx, and right_gene where pair < color_axis_idx
        flow_df$.color <- ifelse(
            flow_df$pair >= color_axis_idx,
            flow_df$left_gene,
            flow_df$right_gene
        )

        n_colors <- length(axis_levels[[color_by]])
        if (is.null(palette)) palette <- .tcrdistR_palette(n_colors)
        names(palette) <- axis_levels[[color_by]]
    }

    # -- build plot ------------------------------------------------------------
    p <- ggplot2::ggplot()

    if (use_fill_aes) {
        p <- p + ggplot2::geom_polygon(
            data = flow_df,
            ggplot2::aes(x = .data$x, y = .data$y,
                         group = .data$.gid,
                         fill = .data$.color),
            alpha = alpha
        ) +
        ggplot2::scale_fill_manual(values = palette, name = color_by)
    } else {
        p <- p + ggplot2::geom_polygon(
            data = flow_df,
            ggplot2::aes(x = .data$x, y = .data$y,
                         group = .data$.gid),
            fill = fill_color, alpha = alpha
        )
    }

    p <- p +
        ggplot2::geom_rect(
            data = strata_df,
            ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                         ymin = .data$ymin, ymax = .data$ymax),
            fill = "grey90", color = "grey40", linewidth = 0.3
        ) +
        ggplot2::geom_text(
            data = strata_df[strata_df$ymax - strata_df$ymin > 0, ],
            ggplot2::aes(
                x = (.data$xmin + .data$xmax) / 2,
                y = (.data$ymin + .data$ymax) / 2,
                label = .data$label
            ),
            size = 2.5, fontface = "plain"
        ) +
        ggplot2::scale_x_continuous(
            breaks = seq_along(axes), labels = axes,
            expand = c(0.08, 0.08)
        ) +
        ggplot2::labs(y = "Count", title = title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
        )

    p
}


# ---------------------------------------------------------------------------
# plot_cdr3_length  (exported)
# ---------------------------------------------------------------------------

#' Plot CDR3 length distribution
#'
#' Draws a bar chart of CDR3 amino-acid sequence lengths for the alpha chain,
#' beta chain, or both (faceted side-by-side).
#'
#' @param tcr_df Data.frame with \code{cdr3a} and/or \code{cdr3b} columns.
#' @param chain Character string. Which chain(s) to plot:
#'   \code{"alpha"}, \code{"beta"}, or \code{"both"} (default). When
#'   \code{"both"}, the plot is faceted by chain.
#' @param max_lengths Integer. Maximum number of distinct lengths to display.
#'   Rare extremes are dropped. Default \code{30L}.
#' @param title Optional character string. Plot title.
#' @param tcr_rep A \code{\linkS4class{TCRrep}} object.  If provided, the
#'   clone data.frame is used as \code{tcr_df}.  Cannot be combined with
#'   \code{tcr_df}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' data(dash)
#' plot_cdr3_length(dash, chain = "both")
#' plot_cdr3_length(dash, chain = "beta", title = "Beta CDR3 lengths")
#' }
#'
#' @seealso \code{\link{plot_gene_usage}}, \code{\link{plot_cdr3_logo}}
#' @export
plot_cdr3_length <- function(tcr_df = NULL, chain = c("alpha", "beta", "both"),
                              max_lengths = 30L, title = NULL,
                              tcr_rep = NULL) {
    .check_ggplot2("plot_cdr3_length()")

    if (!is.null(tcr_rep)) {
        if (!is.null(tcr_df))
            stop("provide 'tcr_rep' or 'tcr_df', not both", call. = FALSE)
        tcr_df <- .extract_from_tcr_rep(tcr_rep)$clone_df
    }

    chain <- match.arg(chain)

    build_lengths <- function(seqs, label) {
        seqs <- seqs[!is.na(seqs) & nzchar(seqs)]
        if (length(seqs) == 0L) return(NULL)
        lens <- nchar(seqs)
        data.frame(length = lens, chain = label, stringsAsFactors = FALSE)
    }

    dfs <- list()
    if (chain %in% c("alpha", "both") && "cdr3a" %in% colnames(tcr_df)) {
        dfs$alpha <- build_lengths(tcr_df$cdr3a, "CDR3\u03b1")
    }
    if (chain %in% c("beta", "both") && "cdr3b" %in% colnames(tcr_df)) {
        dfs$beta <- build_lengths(tcr_df$cdr3b, "CDR3\u03b2")
    }

    df <- do.call(rbind, dfs)
    if (is.null(df) || nrow(df) == 0L) {
        return(ggplot2::ggplot() + ggplot2::theme_void() +
                   ggplot2::ggtitle(title %||% "No CDR3 data"))
    }

    tally <- as.data.frame(table(df$length, df$chain), stringsAsFactors = FALSE)
    colnames(tally) <- c("length", "chain", "count")
    tally$length <- as.integer(tally$length)
    tally <- tally[tally$count > 0L, ]

    if (is.null(title)) {
        title <- "CDR3 length distribution"
    }

    p <- ggplot2::ggplot(tally, ggplot2::aes(x = .data$length, y = .data$count)) +
        ggplot2::geom_col(fill = "#1f77b4", width = 0.7) +
        ggplot2::labs(x = "CDR3 length (aa)", y = "Count", title = title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
        )

    if (chain == "both") {
        p <- p + ggplot2::facet_wrap(~ chain, scales = "free_y")
    }

    p
}


# ---------------------------------------------------------------------------
# .render_gene_text_raster  (internal helper)
# ---------------------------------------------------------------------------

#' Render a gene name label as a cropped RGBA raster
#'
#' Creates a temporary PNG, draws bold text at a large font size, then reads
#' the image back and crops to the bounding box of non-transparent pixels.
#' The resulting array can be passed to \code{ggplot2::annotation_raster()}.
#'
#' @param label Character string. The text to render.
#' @param col Character string. Text colour. Default \code{"black"}.
#' @param width_px Integer. Pixel width of the temporary canvas.
#'   Default \code{2400L}.
#' @param height_px Integer. Pixel height of the temporary canvas.
#'   Default \code{480L}.
#' @return A numeric array (rows x cols x 4 RGBA channels) suitable for
#'   \code{annotation_raster()}.
#' @keywords internal
#' @noRd
.render_gene_text_raster <- function(label, col = "black",
                                      width_px = 2400L, height_px = 480L) {
    if (!requireNamespace("png", quietly = TRUE)) {
        stop("Package 'png' is required for gene logo rendering. ",
             "Install it with: install.packages(\"png\")",
             call. = FALSE)
    }

    tf <- tempfile(fileext = ".png")
    on.exit(unlink(tf), add = TRUE)
    grDevices::png(tf, width = width_px, height = height_px,
                   bg = "transparent")
    grid::grid.newpage()
    grid::grid.text(label,
        gp = grid::gpar(fontface = "bold", col = col, fontsize = 320))
    grDevices::dev.off()
    img <- png::readPNG(tf)

    # Crop to non-transparent content bounding box
    alpha <- img[, , 4]
    rows <- which(rowSums(alpha) > 0)
    cols <- which(colSums(alpha) > 0)
    if (length(rows) == 0L || length(cols) == 0L) return(img)
    img[min(rows):max(rows), min(cols):max(cols), , drop = FALSE]
}


# ---------------------------------------------------------------------------
# plot_vj_gene_logo  (exported)
# ---------------------------------------------------------------------------

#' Plot a V/J gene usage logo
#'
#' Renders gene names as stacked text glyphs with heights proportional to
#' their frequency --- a "sequence logo" style for gene usage. Each gene
#' name is rendered as a coloured raster image and vertically stacked so
#' that more frequent genes occupy more vertical space.
#'
#' @param genes Character vector. V-gene or J-gene allele names
#'   (e.g., \code{"TRAV1-2*01"}).
#' @param organism Character string. Organism identifier
#'   (e.g., \code{"human"}, \code{"mouse"}).
#' @param gene_type Character string. \code{"V"} or \code{"J"}.
#' @param chain Character string. \code{"alpha"} or \code{"beta"}.
#' @param max_genes Integer. Maximum number of genes to display.
#'   Default \code{10L}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' data(dash)
#' pa <- dash[dash$epitope == "PA", ]
#' plot_vj_gene_logo(pa$vb, organism = "mouse", gene_type = "V",
#'                   chain = "beta")
#' plot_vj_gene_logo(pa$ja, organism = "mouse", gene_type = "J",
#'                   chain = "alpha")
#' }
#'
#' @seealso \code{\link{plot_gene_usage}}, \code{\link{plot_tcr_logo_panel}}
#' @export
plot_vj_gene_logo <- function(genes,
                               organism,
                               gene_type = c("V", "J"),
                               chain = c("alpha", "beta"),
                               max_genes = 10L) {
    .check_ggplot2("plot_vj_gene_logo()")

    gene_type <- match.arg(gene_type)
    chain <- match.arg(chain)
    max_genes <- as.integer(max_genes)

    stopifnot(
        is.character(genes),
        length(genes) >= 1L
    )

    # ---- Map alleles to count_rep names --------------------------------------
    organism_genes <- load_gene_database(organism)

    count_reps <- vapply(genes, function(g) {
        entry <- organism_genes[[g]]
        if (!is.null(entry)) {
            entry$count_rep
        } else {
            trim_allele_to_gene(g)
        }
    }, character(1L), USE.NAMES = FALSE)

    # ---- Compute frequencies and proportions ---------------------------------
    freq_table <- sort(table(count_reps), decreasing = TRUE)
    gene_names <- names(freq_table)
    proportions <- as.numeric(freq_table) / length(count_reps)

    # Trim to max_genes
    n_show <- min(length(gene_names), max_genes)
    gene_names  <- gene_names[seq_len(n_show)]
    proportions <- proportions[seq_len(n_show)]

    # ---- Trim display names --------------------------------------------------
    prefixes <- c("TRAV", "TRAJ", "TRBV", "TRBJ",
                   "TRGV", "TRGJ", "TRDV", "TRDJ",
                   "IGHV", "IGHJ", "IGLV", "IGLJ", "IGKV", "IGKJ")
    display_names <- gene_names
    for (pfx in prefixes) {
        mask <- startsWith(display_names, pfx)
        display_names[mask] <- substring(display_names[mask],
                                          nchar(pfx) + 1L)
    }

    # ---- Colour palette ------------------------------------------------------
    pal <- .tcrdistR_palette(n_show)

    # ---- Stack positions (most frequent at bottom) ---------------------------
    y_top <- cumsum(proportions)
    y_bottom <- c(0, y_top[-n_show])

    # ---- Title ---------------------------------------------------------------
    chain_letter <- if (chain == "alpha") "\u03b1" else "\u03b2"
    plot_title <- paste0(gene_type, chain_letter)

    # ---- Build ggplot with raster gene glyphs --------------------------------
    p <- ggplot2::ggplot() +
        ggplot2::scale_x_continuous(limits = c(0, 1),
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits = c(0, max(y_top)),
                                    expand = ggplot2::expansion(
                                        mult = c(0, 0.02))) +
        ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 10,
                                                face = "bold"),
            plot.margin = ggplot2::margin(2, 2, 2, 2)
        )

    for (i in seq_len(n_show)) {
        img <- .render_gene_text_raster(display_names[i], col = pal[i])
        p <- p + ggplot2::annotation_raster(img,
            xmin = 0.05, xmax = 0.95,
            ymin = y_bottom[i], ymax = y_top[i],
            interpolate = TRUE)
    }

    p
}


# ---------------------------------------------------------------------------
# compute_nucseq_src  (exported)
# ---------------------------------------------------------------------------

#' Compute nucleotide source annotations for CDR3 sequences
#'
#' For each TCR, calls \code{.analyze_junction()} to determine the V/N/D/J
#' origin of each nucleotide in the CDR3 region. The result can be passed
#' to \code{\link{plot_cdr3_logo}} via its \code{nucseq_src} parameter to
#' display junction bars showing the rearrangement structure.
#'
#' @param tcrs Data.frame with columns \code{va}, \code{ja}, \code{cdr3a},
#'   \code{cdr3a_nucseq} (for alpha chain) or \code{vb}, \code{jb},
#'   \code{cdr3b}, \code{cdr3b_nucseq} (for beta chain).
#' @param organism Character string. Organism identifier
#'   (e.g., \code{"human"}, \code{"mouse"}).
#' @param chain Character string. \code{"alpha"} or \code{"beta"}.
#'
#' @return A list of character vectors (one per TCR). Each vector has
#'   length \code{nchar(cdr3_nucseq)} with elements from
#'   \code{c("V", "N", "J")} for alpha chain, or
#'   \code{c("V", "N1", "D", "N2", "J")} for beta chain. Returns
#'   \code{NULL} for TCRs where junction analysis fails.
#'
#' @examples
#' \dontrun{
#' data(dash)
#' pa <- dash[dash$epitope == "PA", ][1:20, ]
#' src <- compute_nucseq_src(pa, organism = "mouse", chain = "beta")
#' plot_cdr3_logo(pa$cdr3b, chain = "beta", nucseq_src = src,
#'                show_junction_bars = TRUE)
#' }
#'
#' @seealso \code{\link{plot_cdr3_logo}}, \code{\link{plot_junction_bars}},
#'   \code{\link{plot_tcr_logo_panel}}
#' @export
compute_nucseq_src <- function(tcrs, organism,
                                chain = c("alpha", "beta")) {
    chain <- match.arg(chain)
    n <- nrow(tcrs)
    result <- vector("list", n)

    for (i in seq_len(n)) {
        if (chain == "alpha") {
            v_gene <- tcrs$va[i]
            j_gene <- tcrs$ja[i]
            cdr3   <- tcrs$cdr3a[i]
            nucseq <- tcrs$cdr3a_nucseq[i]
        } else {
            v_gene <- tcrs$vb[i]
            j_gene <- tcrs$jb[i]
            cdr3   <- tcrs$cdr3b[i]
            nucseq <- tcrs$cdr3b_nucseq[i]
        }

        if (is.na(nucseq) || nucseq == "" || is.na(cdr3) || cdr3 == "") {
            next
        }

        res <- tryCatch(
            .analyze_junction(organism, v_gene, j_gene, cdr3,
                              tolower(nucseq)),
            error = function(e) NULL
        )

        if (is.null(res)) next

        src_chars <- strsplit(res$cdr3_nucseq_src, "", fixed = TRUE)[[1L]]

        # For beta chain: convert N -> N1 (before D) or N2 (after D)
        if (chain == "beta") {
            seen_d <- FALSE
            for (j in seq_along(src_chars)) {
                if (src_chars[j] == "D") {
                    seen_d <- TRUE
                } else if (src_chars[j] == "N") {
                    src_chars[j] <- if (seen_d) "N2" else "N1"
                }
            }
        }

        result[[i]] <- src_chars
    }

    result
}


# ---------------------------------------------------------------------------
# plot_tcr_logo_panel  (exported)
# ---------------------------------------------------------------------------

#' Plot a composite TCR rearrangement logo panel
#'
#' Arranges V-gene logos, CDR3 sequence logos (with optional junction bars),
#' and J-gene logos for both alpha and beta chains in a single composite
#' panel. This provides a comprehensive view of TCR rearrangement structure.
#'
#' @param tcrs Data.frame with at least columns \code{va}, \code{ja},
#'   \code{cdr3a}, \code{vb}, \code{jb}, \code{cdr3b}. For junction bars,
#'   also requires \code{cdr3a_nucseq} and \code{cdr3b_nucseq}.
#' @param organism Character string. Organism identifier
#'   (e.g., \code{"human"}, \code{"mouse"}).
#' @param show_junction_bars Logical. If \code{TRUE} (default) and
#'   nucleotide sequence columns are present, display V/N/D/J junction bars
#'   below each CDR3 logo.
#' @param title Optional character string. Overall panel title.
#'
#' @return A \code{patchwork} object.
#'
#' @examples
#' \dontrun{
#' data(dash)
#' pa <- dash[dash$epitope == "PA", ][1:30, ]
#' plot_tcr_logo_panel(pa, organism = "mouse")
#' plot_tcr_logo_panel(pa, organism = "mouse", show_junction_bars = FALSE)
#' }
#'
#' @seealso \code{\link{plot_vj_gene_logo}}, \code{\link{plot_cdr3_logo}},
#'   \code{\link{plot_junction_bars}}, \code{\link{compute_nucseq_src}}
#' @export
plot_tcr_logo_panel <- function(tcrs,
                                 organism,
                                 show_junction_bars = TRUE,
                                 title = NULL) {
    .check_ggplot2("plot_tcr_logo_panel()")
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop("Package 'patchwork' is required for plot_tcr_logo_panel(). ",
             "Install it with: install.packages(\"patchwork\")",
             call. = FALSE)
    }

    stopifnot(is.data.frame(tcrs), nrow(tcrs) >= 1L)

    # ---- Compute junction source annotations if needed -----------------------
    alpha_nucseq_src <- NULL
    beta_nucseq_src <- NULL
    has_nucseq <- show_junction_bars &&
        "cdr3a_nucseq" %in% colnames(tcrs) &&
        "cdr3b_nucseq" %in% colnames(tcrs)
    if (has_nucseq) {
        alpha_nucseq_src <- tryCatch(
            compute_nucseq_src(tcrs, organism, "alpha"),
            error = function(e) NULL
        )
        beta_nucseq_src <- tryCatch(
            compute_nucseq_src(tcrs, organism, "beta"),
            error = function(e) NULL
        )
    }

    panels <- list()
    widths <- numeric(0)

    # ---- Alpha chain: V-gene logo, CDR3 logo, J-gene logo -------------------
    alpha_jbar <- NULL
    alpha_cdr3_idx <- NA_integer_

    panels <- c(panels, list(
        plot_vj_gene_logo(tcrs$va, organism, "V", "alpha")
    ))
    widths <- c(widths, 1.0)

    cdr3a_result <- plot_cdr3_logo(
        tcrs$cdr3a, chain = "alpha",
        nucseq_src = alpha_nucseq_src,
        return_junction_pwm = has_nucseq && !is.null(alpha_nucseq_src)
    )
    if (is.list(cdr3a_result) && !inherits(cdr3a_result, "gg")) {
        panels <- c(panels, list(cdr3a_result$plot))
        if (!is.null(cdr3a_result$junction_pwm)) {
            alpha_jbar <- plot_junction_bars(
                cdr3a_result$junction_pwm, chain = "alpha")
        }
    } else {
        panels <- c(panels, list(cdr3a_result))
    }
    alpha_cdr3_idx <- length(panels)
    widths <- c(widths, 3.0)

    panels <- c(panels, list(
        plot_vj_gene_logo(tcrs$ja, organism, "J", "alpha")
    ))
    widths <- c(widths, 1.0)

    # ---- Beta chain: V-gene logo, CDR3 logo, J-gene logo --------------------
    beta_jbar <- NULL
    beta_cdr3_idx <- NA_integer_

    panels <- c(panels, list(
        plot_vj_gene_logo(tcrs$vb, organism, "V", "beta")
    ))
    widths <- c(widths, 1.0)

    cdr3b_result <- plot_cdr3_logo(
        tcrs$cdr3b, chain = "beta",
        nucseq_src = beta_nucseq_src,
        return_junction_pwm = has_nucseq && !is.null(beta_nucseq_src)
    )
    if (is.list(cdr3b_result) && !inherits(cdr3b_result, "gg")) {
        panels <- c(panels, list(cdr3b_result$plot))
        if (!is.null(cdr3b_result$junction_pwm)) {
            beta_jbar <- plot_junction_bars(
                cdr3b_result$junction_pwm, chain = "beta")
        }
    } else {
        panels <- c(panels, list(cdr3b_result))
    }
    beta_cdr3_idx <- length(panels)
    widths <- c(widths, 3.0)

    panels <- c(panels, list(
        plot_vj_gene_logo(tcrs$jb, organism, "J", "beta")
    ))
    widths <- c(widths, 1.0)

    # ---- Assemble with patchwork ---------------------------------------------
    row1 <- patchwork::wrap_plots(panels, nrow = 1L, widths = widths)

    has_jbars <- !is.null(alpha_jbar) || !is.null(beta_jbar)
    if (has_jbars) {
        # Build junction bar row aligned under CDR3 panels
        spacer <- ggplot2::ggplot() + ggplot2::theme_void()
        jbar_panels <- vector("list", length(panels))
        for (k in seq_along(jbar_panels)) jbar_panels[[k]] <- spacer
        if (!is.null(alpha_jbar)) jbar_panels[[alpha_cdr3_idx]] <- alpha_jbar
        if (!is.null(beta_jbar))  jbar_panels[[beta_cdr3_idx]]  <- beta_jbar

        row2 <- patchwork::wrap_plots(jbar_panels, nrow = 1L, widths = widths)
        result <- patchwork::wrap_plots(row1, row2, ncol = 1L,
                                         heights = c(3, 1))
    } else {
        result <- row1
    }

    if (!is.null(title)) {
        result <- result +
            patchwork::plot_annotation(title = title)
    }

    result
}
