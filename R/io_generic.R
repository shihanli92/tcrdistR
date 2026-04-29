#' @include io_utils.R
NULL

# ---------------------------------------------------------------------------
# read_tcr_table
# ---------------------------------------------------------------------------

#' Read a delimited file containing TCR repertoire data
#'
#' Reads a CSV or TSV file of TCR clonotypes and returns a data.frame with
#' tcrdistR canonical column names (\code{va}, \code{cdr3a}, \code{vb},
#' \code{cdr3b}, and optionally \code{ja}, \code{jb}). Common naming
#' conventions (e.g. tcrdist3/DASH format) are auto-detected when
#' \code{col_map} is not provided.
#'
#' @param file Character string. Path to a delimited text file.
#' @param col_map A named list mapping tcrdistR canonical column names to the
#'   column names used in \code{file}. For example,
#'   \code{list(va = "v_a_gene", cdr3a = "cdr3_a_aa")}. If \code{NULL}
#'   (the default), the function attempts to auto-detect the naming pattern.
#' @param sep Character. Field separator passed to \code{\link[utils]{read.delim}}.
#'   Defaults to \code{"\\t"} (tab-separated).
#' @param normalize_genes Logical. If \code{TRUE} (default), V and J gene
#'   names that lack an allele suffix (e.g. \code{TRAV1-1}) are normalised by
#'   appending \code{*01}.
#' @param ... Additional arguments passed to \code{\link[utils]{read.delim}}.
#'
#' @return A \code{data.frame} with tcrdistR canonical column names. All gene
#'   and CDR3 columns are character vectors (never factors).
#'
#' @examples
#' \donttest{
#' # Read a tcrdist3-style TSV file
#' # df <- read_tcr_table("dash_human.tsv")
#'
#' # Read a CSV with custom column mapping
#' # df <- read_tcr_table("my_data.csv", sep = ",",
#' #     col_map = list(va = "V_alpha", cdr3a = "CDR3_alpha",
#' #                    vb = "V_beta",  cdr3b = "CDR3_beta"))
#' }
#'
#' @seealso \code{\link{read_airr}}, \code{\link{read_adaptive}}, \code{\link{read_10x}}, \code{\link{TCRrep}}
#' @export
read_tcr_table <- function(file, col_map = NULL, sep = "\t",
                           normalize_genes = TRUE, ...) {

    # ---- Validate file path ---------------------------------------------------
    if (!is.character(file) || length(file) != 1L || !nzchar(file)) {
        stop("read_tcr_table: 'file' must be a non-empty character string")
    }
    if (!file.exists(file)) {
        stop(sprintf("read_tcr_table: file not found: %s", file))
    }

    # ---- Read the file --------------------------------------------------------
    df <- utils::read.delim(file, sep = sep, stringsAsFactors = FALSE, ...)

    if (nrow(df) == 0L) {
        stop("read_tcr_table: file contains no data rows")
    }

    # ---- Column mapping -------------------------------------------------------
    if (!is.null(col_map)) {
        if (!is.list(col_map) || is.null(names(col_map))) {
            stop("read_tcr_table: 'col_map' must be a named list")
        }
        df <- .standardize_columns(df, col_map)
    } else {
        df <- .auto_detect_columns(df)
    }

    # ---- Normalize gene names -------------------------------------------------
    if (isTRUE(normalize_genes)) {
        gene_cols <- intersect(c("va", "vb", "ja", "jb"), colnames(df))
        for (col in gene_cols) {
            df[[col]] <- .normalize_gene_names(df[[col]])
        }
    }

    df
}

# ---------------------------------------------------------------------------
# .auto_detect_columns (internal helper)
# ---------------------------------------------------------------------------

#' Auto-detect and rename common TCR column naming patterns
#'
#' @param df A \code{data.frame} as read from a file.
#' @return The \code{data.frame} with columns renamed to canonical names.
#' @keywords internal
.auto_detect_columns <- function(df) {
    cn <- colnames(df)

    # Pattern 1: already canonical
    canonical <- c("va", "cdr3a", "vb", "cdr3b")
    if (all(canonical %in% cn)) {
        return(df)
    }

    # Pattern 2: tcrdist3 / DASH style
    dash_map <- list(
        va    = "v_a_gene",
        cdr3a = "cdr3_a_aa",
        vb    = "v_b_gene",
        cdr3b = "cdr3_b_aa",
        ja    = "j_a_gene",
        jb    = "j_b_gene"
    )
    dash_core <- c("v_a_gene", "cdr3_a_aa", "v_b_gene", "cdr3_b_aa")
    if (all(dash_core %in% cn)) {
        return(.standardize_columns(df, dash_map))
    }

    # No pattern matched
    stop(sprintf(
        paste0("read_tcr_table: could not auto-detect column naming pattern. ",
               "Expected canonical names (va, cdr3a, vb, cdr3b) or ",
               "tcrdist3/DASH names (v_a_gene, cdr3_a_aa, v_b_gene, ",
               "cdr3_b_aa). Found columns: %s. ",
               "Provide an explicit 'col_map' argument."),
        paste(cn, collapse = ", ")
    ))
}
