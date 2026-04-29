# Plot V-gene or J-gene usage frequencies

Horizontal bar chart of gene frequencies from a TCR data.frame column.
Optionally strips allele suffixes (e.g., `"TRAV1-2*01"` becomes
`"TRAV1-2"`).

## Usage

``` r
plot_gene_usage(
  tcr_df,
  gene_col,
  strip_allele = TRUE,
  max_genes = 20L,
  title = NULL
)
```

## Arguments

- tcr_df:

  Data.frame containing TCR data.

- gene_col:

  Character string. Column name to tally (e.g., `"va"`, `"vb"`, `"ja"`).

- strip_allele:

  Logical. If `TRUE` (default), strip allele suffixes (everything after
  `"*"`) before tallying.

- max_genes:

  Integer. Maximum number of genes to display. Default `20L`.

- title:

  Optional plot title.

## Value

A `ggplot` object.

## See also

[`plot_cdr3_logo`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md),
[`plot_tcrdist_dendrogram`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_dendrogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tcr_df <- data.frame(va = c("TRAV1-2*01", "TRAV1-2*02", "TRAV12-1*01",
                             "TRAV1-2*01", "TRAV12-1*01"))
plot_gene_usage(tcr_df, "va", title = "V-alpha usage")
} # }
```
