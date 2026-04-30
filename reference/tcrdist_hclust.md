# Hierarchical clustering of TCRs by TCRdist

Convenience wrapper that computes the pairwise TCRdist matrix and
performs hierarchical clustering via
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

## Usage

``` r
tcrdist_hclust(
  tcr_df = NULL,
  organism = NULL,
  method = "average",
  max_tcrs = 2000L,
  dist_matrix = NULL
)
```

## Arguments

- tcr_df:

  Data.frame with TCR columns (optional if `dist_matrix` is provided).

- organism:

  Character string (`"human"` or `"mouse"`) (optional if `dist_matrix`
  is provided).

- method:

  Clustering method for
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html). Default
  `"average"` (UPGMA).

- max_tcrs:

  Integer. Subsample if N exceeds this. Default `2000L`.

- dist_matrix:

  Optional precomputed distance matrix. If provided, `tcr_df` and
  `organism` are not used for distance computation.

## Value

A named list:

- `hclust`:

  An `hclust` object.

- `dist_matrix`:

  The pairwise distance matrix used.

- `indices`:

  Integer vector of row indices used (after potential subsampling).

## See also

[`cluster_tcrs`](https://shihanli92.github.io/tcrdistR/reference/cluster_tcrs.md),
[`neighborhood_test`](https://shihanli92.github.io/tcrdistR/reference/neighborhood_test.md),
[`plot_tcrdist_dendrogram`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_dendrogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- tcrdist_hclust(tcr_df, "human")
plot(result$hclust)
clusters <- cutree(result$hclust, k = 5)
} # }
```
