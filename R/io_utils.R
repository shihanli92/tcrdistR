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
