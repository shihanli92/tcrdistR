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

    # Pattern 2: scRepertoire (CTgene + CTaa)
    result <- .detect_screpertoire(df, cn)
    if (!is.null(result)) return(result)

    # Pattern 3: scirpy / dandelion (IR_VJ_1_* + IR_VDJ_1_*)
    scirpy_core <- c("IR_VJ_1_v_call", "IR_VDJ_1_v_call")
    if (all(scirpy_core %in% cn)) {
        return(.standardize_columns(df, .scirpy_map()))
    }

    # Pattern 4: tcrdist3 / DASH style
    dash_core <- c("v_a_gene", "cdr3_a_aa", "v_b_gene", "cdr3_b_aa")
    if (all(dash_core %in% cn)) {
        return(.standardize_columns(df, .tcrdist3_map()))
    }

    # No pattern matched
    stop(sprintf(
        paste0("Could not auto-detect column naming pattern. ",
               "Recognized formats: canonical (va, cdr3a, vb, cdr3b), ",
               "scRepertoire (CTgene, CTaa), ",
               "scirpy/dandelion (IR_VJ_1_v_call, IR_VDJ_1_v_call), ",
               "tcrdist3 (v_a_gene, cdr3_a_aa). ",
               "Found columns: %s. ",
               "Provide an explicit 'col_map' argument."),
        paste(cn, collapse = ", ")
    ))
}


# ---------------------------------------------------------------------------
# Format-specific column maps
# ---------------------------------------------------------------------------

#' @keywords internal
.tcrdist3_map <- function() {
    list(
        va    = "v_a_gene",
        cdr3a = "cdr3_a_aa",
        vb    = "v_b_gene",
        cdr3b = "cdr3_b_aa",
        ja    = "j_a_gene",
        jb    = "j_b_gene"
    )
}

#' @keywords internal
.scirpy_map <- function() {
    list(
        va    = "IR_VJ_1_v_call",
        ja    = "IR_VJ_1_j_call",
        cdr3a = "IR_VJ_1_junction_aa",
        vb    = "IR_VDJ_1_v_call",
        jb    = "IR_VDJ_1_j_call",
        cdr3b = "IR_VDJ_1_junction_aa"
    )
}


# ---------------------------------------------------------------------------
# scRepertoire detection and parsing
# ---------------------------------------------------------------------------

#' @keywords internal
.detect_screpertoire <- function(df, cn) {
    # Auto-detect column suffix (_TCR, _BCR, or none)
    suffix <- NULL
    if ("CTgene_TCR" %in% cn) {
        suffix <- "_TCR"
    } else if ("CTgene_BCR" %in% cn) {
        suffix <- "_BCR"
    } else if ("CTgene" %in% cn) {
        suffix <- ""
    }
    if (is.null(suffix)) return(NULL)

    ctgene_col <- paste0("CTgene", suffix)
    ctaa_col   <- paste0("CTaa", suffix)
    ctnt_col   <- paste0("CTnt", suffix)

    if (!ctaa_col %in% cn) return(NULL)

    # Parse gene names and CDR3 sequences
    genes  <- .parse_screpertoire_ctgene(df[[ctgene_col]])
    cdr3aa <- .parse_screpertoire_cdr3(df[[ctaa_col]])

    # Parse nucleotide sequences if available
    has_nt <- ctnt_col %in% cn
    if (has_nt) {
        cdr3nt <- .parse_screpertoire_cdr3(df[[ctnt_col]])
    }

    # Build result: parsed columns + any extra columns from input
    parsed_cols <- c(ctgene_col, ctaa_col)
    if (has_nt) parsed_cols <- c(parsed_cols, ctnt_col)
    extra_cols <- setdiff(cn, parsed_cols)

    result <- data.frame(
        va    = genes$va,
        ja    = genes$ja,
        cdr3a = cdr3aa$alpha,
        vb    = genes$vb,
        jb    = genes$jb,
        cdr3b = cdr3aa$beta,
        stringsAsFactors = FALSE
    )

    if (has_nt) {
        result$cdr3a_nucseq <- cdr3nt$alpha
        result$cdr3b_nucseq <- cdr3nt$beta
    }

    # Carry over extra columns
    if (length(extra_cols) > 0L) {
        result <- cbind(result, df[, extra_cols, drop = FALSE])
    }

    result
}


# ---------------------------------------------------------------------------
# as_tcr_df
# ---------------------------------------------------------------------------

#' Standardize a TCR data.frame to tcrdistR column names
#'
#' Converts a data.frame from common TCR analysis tools to tcrdistR's canonical
#' column names (\code{va}, \code{cdr3a}, \code{vb}, \code{cdr3b}).  Supports
#' auto-detection of column naming patterns from scRepertoire, scirpy,
#' dandelion, and tcrdist3, as well as custom column mappings.
#'
#' @param df A \code{data.frame} containing TCR data.
#' @param col_map Named character vector mapping tcrdistR canonical names to
#'   the column names in \code{df}.  For example,
#'   \code{c(va = "alpha_v_gene", cdr3a = "alpha_cdr3")}.
#'   If \code{NULL}, the format is auto-detected.
#' @param format Character string or \code{NULL}.  Force a specific format
#'   instead of auto-detecting:
#'   \describe{
#'     \item{\code{"screpertoire"}}{scRepertoire \code{CTgene}/\code{CTaa}
#'       concatenated columns.}
#'     \item{\code{"scirpy"}, \code{"dandelion"}}{scirpy/dandelion
#'       \code{IR_VJ_1_*}/\code{IR_VDJ_1_*} columns.}
#'     \item{\code{"tcrdist3"}}{tcrdist3 \code{v_a_gene}/\code{cdr3_a_aa}
#'       columns.}
#'   }
#'   If \code{NULL} (default), auto-detection is used.
#' @param normalize_genes Logical.  If \code{TRUE} (default), gene names
#'   without an allele suffix get \code{*01} appended.
#' @param drop_incomplete Logical.  If \code{TRUE} (default), rows with
#'   \code{NA} or empty strings in required chain columns (\code{va},
#'   \code{cdr3a}, \code{vb}, \code{cdr3b}) are dropped with a message.
#'
#' @return A \code{data.frame} with tcrdistR canonical column names.  Extra
#'   columns from the input are preserved.
#'
#' @examples
#' # Custom column mapping
#' df <- data.frame(
#'   alpha_v = "TRAV1-1", alpha_cdr3 = "CAVRDSSYKLIF",
#'   beta_v  = "TRBV5-1", beta_cdr3  = "CASSIRSSYEQYF"
#' )
#' as_tcr_df(df, col_map = c(va = "alpha_v", cdr3a = "alpha_cdr3",
#'                            vb = "beta_v",  cdr3b = "beta_cdr3"))
#'
#' # scirpy / dandelion format (auto-detected)
#' df <- data.frame(
#'   IR_VJ_1_v_call = "TRAV1-1", IR_VJ_1_junction_aa = "CAVRDSSYKLIF",
#'   IR_VDJ_1_v_call = "TRBV5-1", IR_VDJ_1_junction_aa = "CASSIRSSYEQYF"
#' )
#' as_tcr_df(df)
#'
#' # scRepertoire format (auto-detected)
#' df <- data.frame(
#'   CTgene = "TRAV1-1.TRAJ33.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
#'   CTaa   = "CAVRDSSYKLIF_CASSIRSSYEQYF"
#' )
#' as_tcr_df(df)
#'
#' @seealso \code{\link{read_tcr_table}}, \code{\link{TCRrep}}
#' @export
as_tcr_df <- function(df,
                       col_map = NULL,
                       format = NULL,
                       normalize_genes = TRUE,
                       drop_incomplete = TRUE) {
    if (!is.data.frame(df)) {
        stop("as_tcr_df: 'df' must be a data.frame", call. = FALSE)
    }
    if (nrow(df) == 0L) {
        stop("as_tcr_df: 'df' has zero rows", call. = FALSE)
    }

    # ---- Resolution: col_map > format > auto-detect --------------------------
    if (!is.null(col_map)) {
        # Custom mapping
        if (is.null(names(col_map)) || !is.character(col_map)) {
            stop("as_tcr_df: 'col_map' must be a named character vector, ",
                 "e.g. c(va = \"v_alpha\", cdr3a = \"cdr3_alpha\")",
                 call. = FALSE)
        }
        # Convert to list for .standardize_columns()
        df <- .standardize_columns(df, as.list(col_map))
        fmt_name <- "custom"
    } else if (!is.null(format)) {
        format <- match.arg(format, c("screpertoire", "scirpy", "dandelion",
                                       "tcrdist3"))
        if (format == "screpertoire") {
            result <- .detect_screpertoire(df, colnames(df))
            if (is.null(result)) {
                stop("as_tcr_df: format='screpertoire' but CTgene/CTaa ",
                     "columns not found", call. = FALSE)
            }
            df <- result
        } else if (format %in% c("scirpy", "dandelion")) {
            df <- .standardize_columns(df, .scirpy_map())
        } else if (format == "tcrdist3") {
            df <- .standardize_columns(df, .tcrdist3_map())
        }
        fmt_name <- format
    } else {
        # Auto-detect
        fmt_name <- .detect_format_name(df)
        df <- .auto_detect_columns(df)
    }

    message("as_tcr_df: detected format '", fmt_name, "'")

    # ---- Normalize gene names ------------------------------------------------
    if (isTRUE(normalize_genes)) {
        gene_cols <- intersect(c("va", "vb", "ja", "jb"), colnames(df))
        for (col in gene_cols) {
            df[[col]] <- .normalize_gene_names(df[[col]])
        }
    }

    # ---- Drop incomplete rows ------------------------------------------------
    if (isTRUE(drop_incomplete)) {
        required <- intersect(c("va", "cdr3a", "vb", "cdr3b"), colnames(df))
        if (length(required) > 0L) {
            complete <- rowSums(
                is.na(df[, required, drop = FALSE]) |
                    df[, required, drop = FALSE] == ""
            ) == 0L
            n_dropped <- sum(!complete)
            if (n_dropped > 0L) {
                message("as_tcr_df: dropped ", n_dropped,
                        " rows with incomplete TCR data (",
                        sum(complete), " retained)")
                df <- df[complete, , drop = FALSE]
            }
        }
    }

    # Coerce factors to character in chain columns
    chain_cols <- intersect(c("va", "ja", "cdr3a", "vb", "jb", "cdr3b"),
                            colnames(df))
    for (col in chain_cols) {
        if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
    }

    rownames(df) <- NULL
    df
}


#' Identify the format name without transforming the data.frame
#' @keywords internal
.detect_format_name <- function(df) {
    cn <- colnames(df)
    canonical <- c("va", "cdr3a", "vb", "cdr3b")
    if (all(canonical %in% cn)) return("canonical")

    if (any(c("CTgene", "CTgene_TCR", "CTgene_BCR") %in% cn)) {
        return("screpertoire")
    }

    scirpy_core <- c("IR_VJ_1_v_call", "IR_VDJ_1_v_call")
    if (all(scirpy_core %in% cn)) return("scirpy")

    dash_core <- c("v_a_gene", "cdr3_a_aa", "v_b_gene", "cdr3_b_aa")
    if (all(dash_core %in% cn)) return("tcrdist3")

    "unknown"
}
