# Summarize a meta-clonotype's composition

Reports the V/J gene usage, CDR3 length distribution, and subject
representation within a meta-clonotype neighborhood.

## Usage

``` r
summarize_meta_clonotype(tcr_df, center_idx, neighbor_indices)
```

## Arguments

- tcr_df:

  Data.frame with TCR columns.

- center_idx:

  Integer. Row index of the center TCR.

- neighbor_indices:

  Integer vector. Row indices of neighbors.

## Value

A named list:

- `center`:

  Named list of center TCR fields.

- `n_members`:

  Number of members.

- `va_usage`:

  Named integer vector of V-alpha gene counts.

- `vb_usage`:

  Named integer vector of V-beta gene counts.

- `cdr3a_lengths`:

  Integer vector of CDR3-alpha lengths.

- `cdr3b_lengths`:

  Integer vector of CDR3-beta lengths.

## See also

[`find_meta_clonotypes`](https://shihanli92.github.io/tcrdistR/reference/find_meta_clonotypes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
meta <- find_meta_clonotypes(tcr_df, "mouse", radius = 30,
                              subject_col = "subject")
if (nrow(meta) > 0) {
    indices <- as.integer(strsplit(meta$neighbor_indices[1], ",")[[1]])
    summary <- summarize_meta_clonotype(
        tcr_df, meta$center_index[1], indices)
    summary$n_members
    summary$va_usage
}
} # }
```
