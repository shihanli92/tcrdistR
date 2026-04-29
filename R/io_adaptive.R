# Reader for Adaptive Biotechnologies ImmunoSeq TSV export files.

# ---------------------------------------------------------------------------
# .convert_adaptive_gene
# ---------------------------------------------------------------------------

#' Convert Adaptive Biotechnologies gene names to IMGT format
#'
#' Transforms gene names from Adaptive's naming convention (e.g.
#' \code{TCRBV05-01*01}) to standard IMGT format (e.g. \code{TRBV5-1*01}).
#' Handles the prefix change (\code{TCR} to \code{TR}), and leading-zero
#' stripping from both family and sub-family numbers.
#'
#' @param gene_names Character vector of Adaptive-format gene names.
#' @return Character vector of the same length with IMGT-formatted names.
#'   NAs and unresolved values are returned as \code{NA_character_}.
#' @keywords internal
.convert_adaptive_gene <- function(gene_names) {
    if (!is.character(gene_names)) {
        stop(".convert_adaptive_gene: 'gene_names' must be a character vector")
    }

    result <- character(length(gene_names))

    for (i in seq_along(gene_names)) {
        g <- gene_names[i]

        # Handle NA, empty, and unresolved
        if (is.na(g) || !nzchar(g) || tolower(g) == "unresolved") {
            result[i] <- NA_character_
            next
        }

        # Step 1: Replace TCR[AB][VDJ] prefix with TR[AB][VDJ]
        g <- sub("^TCR", "TR", g)

        # Step 2: Strip leading zeros from family number (e.g. TRBV05 -> TRBV5)
        g <- sub("^(TR[AB][VDJ])0*(\\d+)", "\\1\\2", g)

        # Step 3: Strip leading zeros from sub-family number (e.g. -01 -> -1)
        g <- sub("-0*(\\d+)", "-\\1", g)

        result[i] <- g
    }

    result
}


# ---------------------------------------------------------------------------
# read_adaptive
# ---------------------------------------------------------------------------

#' Read Adaptive Biotechnologies ImmunoSeq TSV file
#'
#' Imports a beta-chain TCR repertoire from an Adaptive ImmunoSeq export
#' file (TSV). Gene names are converted from Adaptive naming to IMGT format
#' and optionally normalized to include an allele suffix.
#'
#' Two column naming conventions are supported automatically:
#' \describe{
#'   \item{Convention 1 (newer)}{\code{vGeneName}, \code{jGeneName},
#'     \code{aminoAcid}, \code{count}}
#'   \item{Convention 2 (older)}{\code{v_resolved} (or \code{vFamilyName}),
#'     \code{j_resolved} (or \code{jFamilyName}), \code{aminoAcid},
#'     \code{count} / \code{templates} / \code{reads}}
#' }
#'
#' @param file Character string. Path to an Adaptive ImmunoSeq TSV file.
#' @param normalize_genes Logical. If \code{TRUE} (the default), gene names
#'   that lack an allele suffix (\code{*XX}) are appended with \code{*01}
#'   after Adaptive-to-IMGT conversion.
#' @return A \code{data.frame} with columns \code{vb}, \code{jb},
#'   \code{cdr3b}, and optionally \code{count}. Gene names are in IMGT
#'   format. Rows with empty or NA CDR3 sequences are removed.
#' @examples
#' \donttest{
#' # df <- read_adaptive("sample_immunoseq.tsv")
#' # head(df)
#' }
#' @export
read_adaptive <- function(file, normalize_genes = TRUE) {
    # ---- Validate inputs ------------------------------------------------------
    if (!is.character(file) || length(file) != 1L || !nzchar(file)) {
        stop("read_adaptive: 'file' must be a non-empty character string")
    }
    if (!file.exists(file)) {
        stop(sprintf("read_adaptive: file not found: %s", file))
    }

    # ---- Read TSV -------------------------------------------------------------
    raw <- utils::read.delim(file, stringsAsFactors = FALSE)
    cols <- colnames(raw)

    # ---- Detect V/J column convention -----------------------------------------
    if ("vGeneName" %in% cols) {
        v_col <- "vGeneName"
        j_col <- "jGeneName"
    } else if ("v_resolved" %in% cols) {
        v_col <- "v_resolved"
        j_col <- "j_resolved"
    } else if ("vFamilyName" %in% cols) {
        v_col <- "vFamilyName"
        j_col <- "jFamilyName"
    } else {
        stop(sprintf(
            "read_adaptive: could not detect column convention. Available columns: %s",
            paste(cols, collapse = ", ")
        ))
    }

    # ---- Detect count column (optional) ---------------------------------------
    count_col <- NULL
    for (candidate in c("count", "templates", "reads")) {
        if (candidate %in% cols) {
            count_col <- candidate
            break
        }
    }

    # ---- Build output data.frame ----------------------------------------------
    out <- data.frame(
        vb    = raw[[v_col]],
        jb    = raw[[j_col]],
        cdr3b = raw[["aminoAcid"]],
        stringsAsFactors = FALSE
    )
    if (!is.null(count_col)) {
        out$count <- raw[[count_col]]
    }

    # ---- Filter rows with missing CDR3 ---------------------------------------
    keep <- !is.na(out$cdr3b) & nzchar(out$cdr3b)
    out <- out[keep, , drop = FALSE]
    rownames(out) <- NULL

    # ---- Convert gene names ---------------------------------------------------
    out$vb <- .convert_adaptive_gene(out$vb)
    out$jb <- .convert_adaptive_gene(out$jb)

    if (normalize_genes) {
        out$vb <- .normalize_gene_names(out$vb)
        out$jb <- .normalize_gene_names(out$jb)
    }

    out
}
