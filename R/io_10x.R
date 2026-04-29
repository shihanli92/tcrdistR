# 10X Genomics VDJ format reader

#' Read 10X Genomics VDJ contig annotations
#'
#' Reads a 10X Genomics cellranger VDJ output file (typically
#' \code{filtered_contig_annotations.csv}) and returns a paired alpha-beta TCR
#' data frame suitable for \code{\link{TCRrep}}. Each row in the 10X format
#' represents a single contig; this function filters, deduplicates, and pairs
#' alpha and beta chains by cell barcode.
#'
#' When a cell has multiple contigs for the same chain, the contig with the
#' highest UMI count is retained (if the \code{umis} column is present);
#' otherwise the first contig is kept.
#'
#' @param file Character string. Path to a 10X VDJ CSV file (e.g.
#'   \code{filtered_contig_annotations.csv}).
#' @param pair_by Character string. Column name used to pair TRA and TRB
#'   contigs from the same cell. Defaults to \code{"barcode"}.
#' @param productive_only Logical. If \code{TRUE} (default), only productive
#'   contigs are retained (\code{productive == "True"}).
#' @param normalize_genes Logical. If \code{TRUE} (default), gene names lacking
#'   an allele suffix (e.g. \code{"TRAV1-1"}) are appended with \code{*01}.
#'
#' @return A \code{data.frame} with columns \code{barcode}, \code{va},
#'   \code{ja}, \code{cdr3a}, \code{vb}, \code{jb}, \code{cdr3b}. Only cells
#'   with both a TRA and TRB contig are included (inner join).
#'
#' @examples
#' \donttest{
#' # df <- read_10x("filtered_contig_annotations.csv")
#' # obj <- TCRrep(df, organism = "human")
#' }
#'
#' @export
read_10x <- function(file,
                     pair_by         = "barcode",
                     productive_only = TRUE,
                     normalize_genes = TRUE) {

    # ---- Validate file exists -----------------------------------------------
    if (!is.character(file) || length(file) != 1L || !nzchar(file)) {
        stop("read_10x: 'file' must be a non-empty character string")
    }
    if (!file.exists(file)) {
        stop(sprintf("read_10x: file not found: %s", file))
    }

    # ---- Read CSV -----------------------------------------------------------
    df <- utils::read.csv(file, stringsAsFactors = FALSE)

    # ---- Validate required columns ------------------------------------------
    required <- c(pair_by, "chain", "v_gene", "j_gene", "cdr3", "productive")
    missing_cols <- setdiff(required, colnames(df))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "read_10x: missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    # ---- Filter: is_cell, full_length, productive ---------------------------
    if ("is_cell" %in% colnames(df)) {
        df <- df[df$is_cell == "True", , drop = FALSE]
    }
    if ("full_length" %in% colnames(df)) {
        df <- df[df$full_length == "True", , drop = FALSE]
    }
    if (isTRUE(productive_only)) {
        df <- df[df$productive == "True", , drop = FALSE]
    }

    if (nrow(df) == 0L) {
        stop("read_10x: no rows remain after filtering")
    }

    # ---- Split by chain -----------------------------------------------------
    tra <- df[df$chain == "TRA", , drop = FALSE]
    trb <- df[df$chain == "TRB", , drop = FALSE]

    if (nrow(tra) == 0L) stop("read_10x: no TRA contigs found")
    if (nrow(trb) == 0L) stop("read_10x: no TRB contigs found")

    # ---- Deduplicate: one contig per cell per chain -------------------------
    .dedup_chain <- function(chain_df, id_col) {
        has_umis <- "umis" %in% colnames(chain_df)
        if (has_umis) {
            chain_df <- chain_df[order(chain_df[[id_col]],
                                       -chain_df$umis), , drop = FALSE]
        }
        chain_df[!duplicated(chain_df[[id_col]]), , drop = FALSE]
    }

    tra <- .dedup_chain(tra, pair_by)
    trb <- .dedup_chain(trb, pair_by)

    # ---- Rename and subset columns ------------------------------------------
    tra_out <- data.frame(
        barcode = tra[[pair_by]],
        va      = tra$v_gene,
        ja      = tra$j_gene,
        cdr3a   = tra$cdr3,
        stringsAsFactors = FALSE
    )
    trb_out <- data.frame(
        barcode = trb[[pair_by]],
        vb      = trb$v_gene,
        jb      = trb$j_gene,
        cdr3b   = trb$cdr3,
        stringsAsFactors = FALSE
    )

    # ---- Inner join by barcode ----------------------------------------------
    paired <- merge(tra_out, trb_out, by = "barcode")

    if (nrow(paired) == 0L) {
        stop("read_10x: no paired TRA/TRB contigs found after joining")
    }

    # ---- Normalize gene names -----------------------------------------------
    if (isTRUE(normalize_genes)) {
        gene_cols <- c("va", "ja", "vb", "jb")
        for (col in gene_cols) {
            paired[[col]] <- .normalize_gene_names(paired[[col]])
        }
    }

    # ---- Return standardized column order -----------------------------------
    paired[, c("barcode", "va", "ja", "cdr3a", "vb", "jb", "cdr3b"),
           drop = FALSE]
}
