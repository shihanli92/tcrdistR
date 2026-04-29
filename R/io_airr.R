# AIRR Community Standard format reader

#' Read AIRR Rearrangement TSV format
#'
#' Reads a file in the AIRR Community Standard Rearrangement TSV format and
#' returns a paired alpha-beta TCR data frame suitable for \code{\link{TCRrep}}.
#' Each row in the AIRR format represents a single rearrangement; this function
#' pairs alpha and beta chains by a shared identifier (typically \code{cell_id}).
#'
#' Multi-value gene calls (comma-separated, e.g. \code{"TRAV1-1*01,TRAV1-2*01"})
#' are resolved by taking the first listed allele.
#'
#' @param file Character string. Path to an AIRR Rearrangement TSV file.
#' @param pair_by Character string. Column name used to pair TRA and TRB
#'   rearrangements from the same cell. Defaults to \code{"cell_id"}.
#' @param productive_only Logical. If \code{TRUE} (default), only productive
#'   rearrangements are retained. The \code{productive} column may contain
#'   logical \code{TRUE} or character \code{"T"}/\code{"TRUE"}.
#' @param normalize_genes Logical. If \code{TRUE} (default), gene names lacking
#'   an allele suffix (e.g. \code{"TRAV1-1"}) are appended with \code{*01}.
#'
#' @return A \code{data.frame} with columns \code{cell_id}, \code{va},
#'   \code{ja}, \code{cdr3a}, \code{vb}, \code{jb}, \code{cdr3b}. Only cells
#'   with both a TRA and TRB rearrangement are included (inner join).
#'
#' @examples
#' \donttest{
#' # df <- read_airr("rearrangements.tsv")
#' # obj <- TCRrep(df, organism = "human")
#' }
#'
#' @export
read_airr <- function(file,
                      pair_by         = "cell_id",
                      productive_only = TRUE,
                      normalize_genes = TRUE) {

    # ---- Validate file exists -----------------------------------------------
    if (!is.character(file) || length(file) != 1L || !nzchar(file)) {
        stop("read_airr: 'file' must be a non-empty character string")
    }
    if (!file.exists(file)) {
        stop(sprintf("read_airr: file not found: %s", file))
    }

    # ---- Read TSV -----------------------------------------------------------
    df <- utils::read.delim(file, stringsAsFactors = FALSE)

    # ---- Validate required columns ------------------------------------------
    required <- c("v_call", "j_call", "junction_aa", "locus", pair_by)
    missing_cols <- setdiff(required, colnames(df))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "read_airr: missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    # ---- Filter to productive rearrangements --------------------------------
    if (isTRUE(productive_only)) {
        if ("productive" %in% colnames(df)) {
            keep <- df$productive %in% c(TRUE, "T", "TRUE", "true", "True")
            df <- df[keep, , drop = FALSE]
        }
    }

    if (nrow(df) == 0L) {
        stop("read_airr: no rows remain after filtering")
    }

    # ---- Split by locus -----------------------------------------------------
    tra <- df[df$locus == "TRA", , drop = FALSE]
    trb <- df[df$locus == "TRB", , drop = FALSE]

    if (nrow(tra) == 0L) stop("read_airr: no TRA rearrangements found")
    if (nrow(trb) == 0L) stop("read_airr: no TRB rearrangements found")

    # ---- Resolve multi-value gene calls (take first) ------------------------
    .take_first <- function(x) {
        sub(",.*", "", x)
    }
    tra$v_call <- .take_first(tra$v_call)
    tra$j_call <- .take_first(tra$j_call)
    trb$v_call <- .take_first(trb$v_call)
    trb$j_call <- .take_first(trb$j_call)

    # ---- Rename and subset columns ------------------------------------------
    tra_out <- data.frame(
        cell_id = tra[[pair_by]],
        va      = tra$v_call,
        ja      = tra$j_call,
        cdr3a   = tra$junction_aa,
        stringsAsFactors = FALSE
    )
    trb_out <- data.frame(
        cell_id = trb[[pair_by]],
        vb      = trb$v_call,
        jb      = trb$j_call,
        cdr3b   = trb$junction_aa,
        stringsAsFactors = FALSE
    )

    # ---- Inner join by pair_by ----------------------------------------------
    paired <- merge(tra_out, trb_out, by = "cell_id")

    if (nrow(paired) == 0L) {
        stop("read_airr: no paired TRA/TRB rearrangements found after joining")
    }

    # ---- Normalize gene names -----------------------------------------------
    if (isTRUE(normalize_genes)) {
        gene_cols <- c("va", "ja", "vb", "jb")
        for (col in gene_cols) {
            paired[[col]] <- .normalize_gene_names(paired[[col]])
        }
    }

    # ---- Return standardized column order -----------------------------------
    paired[, c("cell_id", "va", "ja", "cdr3a", "vb", "jb", "cdr3b"),
           drop = FALSE]
}
