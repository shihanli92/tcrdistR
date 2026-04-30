# Plot CDR3 length distribution

Draws a bar chart of CDR3 amino-acid sequence lengths for the alpha
chain, beta chain, or both (faceted side-by-side).

## Usage

``` r
plot_cdr3_length(
  tcr_df,
  chain = c("alpha", "beta", "both"),
  max_lengths = 30L,
  title = NULL
)
```

## Arguments

- tcr_df:

  Data.frame with `cdr3a` and/or `cdr3b` columns.

- chain:

  Character string. Which chain(s) to plot: `"alpha"`, `"beta"`, or
  `"both"` (default). When `"both"`, the plot is faceted by chain.

- max_lengths:

  Integer. Maximum number of distinct lengths to display. Rare extremes
  are dropped. Default `30L`.

- title:

  Optional character string. Plot title.

## Value

A `ggplot` object.

## See also

[`plot_gene_usage`](https://shihanli92.github.io/tcrdistR/reference/plot_gene_usage.md),
[`plot_cdr3_logo`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(dash)
plot_cdr3_length(dash, chain = "both")
plot_cdr3_length(dash, chain = "beta", title = "Beta CDR3 lengths")
} # }
```
