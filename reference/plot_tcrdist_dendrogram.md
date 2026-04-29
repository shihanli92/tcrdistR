# Plot a TCRdist-based dendrogram

Computes pairwise TCRdist distances, performs hierarchical clustering,
and renders a dendrogram with optional leaf coloring.

## Usage

``` r
plot_tcrdist_dendrogram(
  tcr_df,
  organism,
  color_by = NULL,
  method = "average",
  max_tcrs = 500L,
  title = NULL,
  point_size = 2
)
```

## Arguments

- tcr_df:

  Data.frame with TCR columns (`va`, `vb`, `cdr3a`, `cdr3b`).

- organism:

  Character string (`"human"` or `"mouse"`).

- color_by:

  Optional vector for leaf coloring. Numeric gives a viridis gradient;
  character/factor gives a qualitative palette. `NULL` colors all leaves
  gray.

- method:

  Clustering method for
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html). Default
  `"average"` (UPGMA).

- max_tcrs:

  Integer. If `nrow(tcr_df)` exceeds this, subsample. Default `500L`.

- title:

  Optional plot title.

- point_size:

  Numeric. Leaf point size. Default `2`.

## Value

A `ggplot` object.

## See also

[`plot_tcrdist_heatmap`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_heatmap.md),
[`tcrdist_hclust`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_hclust.md),
[`cluster_tcrs`](https://shihanli92.github.io/tcrdistR/reference/cluster_tcrs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot_tcrdist_dendrogram(tcr_df, "human", color_by = tcr_df$cluster)
} # }
```
