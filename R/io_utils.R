# Internal utility functions for TCR repertoire I/O parsers.

# ---------------------------------------------------------------------------
# .normalize_gene_names
# ---------------------------------------------------------------------------

#' Normalize V/J gene names by appending default allele suffix
#'
#' If a gene name lacks a \code{*} allele suffix, appends \code{*01}.
#' NAs and empty strings are passed through unchanged.
#'
#' @param gene_names Character vector of gene names.
#' @return Character vector of the same length with allele suffixes ensured.
#' @keywords internal
.normalize_gene_names <- function(gene_names) {
    if (!is.character(gene_names)) {
        stop(".normalize_gene_names: 'gene_names' must be a character vector")
    }
    needs_allele <- !is.na(gene_names) & nzchar(gene_names) &
        !grepl("*", gene_names, fixed = TRUE)
    gene_names[needs_allele] <- paste0(gene_names[needs_allele], "*01")
    gene_names
}

# ---------------------------------------------------------------------------
# .validate_tcr_df
# ---------------------------------------------------------------------------

#' Validate and clean a TCR data.frame
#'
#' Checks that required columns exist for the specified chain type, ensures
#' no NAs in critical columns, coerces factors to character, and returns
#' the cleaned data.frame.
#'
#' @param df A \code{data.frame} of TCR clonotypes.
#' @param chains Character string. One of \code{"AB"}, \code{"A"}, or
#'   \code{"B"}.
#' @return The cleaned \code{data.frame}.
#' @keywords internal
.validate_tcr_df <- function(df, chains = "AB") {
    if (!is.data.frame(df)) {
        stop(".validate_tcr_df: 'df' must be a data.frame")
    }

    valid_chains <- c("AB", "A", "B")
    if (!chains %in% valid_chains) {
        stop(sprintf(
            ".validate_tcr_df: 'chains' must be one of: %s",
            paste(valid_chains, collapse = ", ")
        ))
    }

    required_cols <- switch(
        chains,
        "AB" = c("va", "cdr3a", "vb", "cdr3b"),
        "A"  = c("va", "cdr3a"),
        "B"  = c("vb", "cdr3b")
    )

    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            ".validate_tcr_df: missing required columns for chains='%s': %s",
            chains, paste(missing_cols, collapse = ", ")
        ))
    }

    # Coerce factors to character
    for (col in required_cols) {
        if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
    }

    # Check for NAs in critical columns
    for (col in required_cols) {
        if (anyNA(df[[col]])) {
            stop(sprintf(
                ".validate_tcr_df: column '%s' contains NA values", col
            ))
        }
    }

    df
}

# ---------------------------------------------------------------------------
# .standardize_columns
# ---------------------------------------------------------------------------

#' Rename columns from source format to tcrdistR canonical names
#'
#' @param df A \code{data.frame}.
#' @param col_map A named list where names are target (canonical) column names
#'   and values are source column names present in \code{df}.
#' @return The \code{data.frame} with matching columns renamed.
#' @keywords internal
.standardize_columns <- function(df, col_map) {
    if (!is.data.frame(df)) {
        stop(".standardize_columns: 'df' must be a data.frame")
    }
    if (!is.list(col_map) || is.null(names(col_map))) {
        stop(".standardize_columns: 'col_map' must be a named list")
    }

    current_names <- colnames(df)
    for (target in names(col_map)) {
        source <- col_map[[target]]
        idx <- match(source, current_names)
        if (!is.na(idx)) {
            current_names[idx] <- target
        }
    }
    colnames(df) <- current_names
    df
}


# ---------------------------------------------------------------------------
# scRepertoire parsing helpers (ported from rconga/R/preprocess.R)
# ---------------------------------------------------------------------------

#' Parse scRepertoire CTgene column into V/J gene columns
#'
#' Parses scRepertoire's concatenated gene format:
#' \code{TRAV.TRAJ.TRAC_TRBV.TRBD.TRBJ.TRBC} (alpha V.J.C underscore
#' beta V.D.J.C).  Multi-chain cells (semicolon-separated) use the first
#' chain per locus.
#'
#' @param ctgene_col Character vector.  The \code{CTgene} column values.
#' @return A \code{data.frame} with columns \code{va}, \code{ja}, \code{vb},
#'   \code{jb}.
#' @keywords internal
.parse_screpertoire_ctgene <- function(ctgene_col) {
    n  <- length(ctgene_col)
    va <- rep(NA_character_, n)
    ja <- rep(NA_character_, n)
    vb <- rep(NA_character_, n)
    jb <- rep(NA_character_, n)

    for (i in seq_len(n)) {
        ct <- ctgene_col[i]
        if (is.na(ct) || ct == "" || ct == "NA") next

        # Split alpha_beta on underscore
        chains <- strsplit(ct, "_", fixed = TRUE)[[1L]]
        if (length(chains) < 2L) next

        # Alpha chain: take first if multi-chain (semicolon-separated)
        alpha_str   <- strsplit(chains[1L], ";", fixed = TRUE)[[1L]][1L]
        alpha_parts <- strsplit(alpha_str, ".", fixed = TRUE)[[1L]]
        # Alpha format: V.J.C (3 parts)
        if (length(alpha_parts) >= 2L) {
            va[i] <- alpha_parts[1L]
            ja[i] <- alpha_parts[2L]
        }

        # Beta chain: take first if multi-chain
        beta_str   <- strsplit(chains[2L], ";", fixed = TRUE)[[1L]][1L]
        beta_parts <- strsplit(beta_str, ".", fixed = TRUE)[[1L]]
        # Beta format: V.D.J.C (4 parts) or V.J.C (3 parts)
        if (length(beta_parts) >= 4L) {
            vb[i] <- beta_parts[1L]
            jb[i] <- beta_parts[3L]
        } else if (length(beta_parts) >= 2L) {
            vb[i] <- beta_parts[1L]
            jb[i] <- beta_parts[2L]
        }
    }

    data.frame(va = va, ja = ja, vb = vb, jb = jb,
               stringsAsFactors = FALSE)
}


#' Parse scRepertoire CTaa/CTnt column into alpha/beta CDR3 columns
#'
#' Parses scRepertoire's concatenated CDR3 format:
#' \code{cdr3a_cdr3b} (underscore-separated).  Multi-chain cells
#' (semicolon-separated) use the first chain per locus.
#'
#' @param cdr3_col Character vector.  The \code{CTaa} or \code{CTnt} column.
#' @return A \code{data.frame} with columns \code{alpha} and \code{beta}.
#' @keywords internal
.parse_screpertoire_cdr3 <- function(cdr3_col) {
    n     <- length(cdr3_col)
    alpha <- rep(NA_character_, n)
    beta  <- rep(NA_character_, n)

    for (i in seq_len(n)) {
        val <- cdr3_col[i]
        if (is.na(val) || val == "" || val == "NA") next

        parts <- strsplit(val, "_", fixed = TRUE)[[1L]]
        if (length(parts) >= 1L) {
            alpha[i] <- strsplit(parts[1L], ";", fixed = TRUE)[[1L]][1L]
        }
        if (length(parts) >= 2L) {
            beta[i] <- strsplit(parts[2L], ";", fixed = TRUE)[[1L]][1L]
        }
    }

    data.frame(alpha = alpha, beta = beta, stringsAsFactors = FALSE)
}
