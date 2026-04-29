# Cluster TCRs by TCRdist-based hierarchical clustering

Computes TCRdist distances, performs hierarchical clustering, and cuts
the tree into groups using either a fixed number of clusters (`k`) or a
height threshold (`h`).

## Usage

``` r
cluster_tcrs(tcr_df, organism, k = NULL, h = NULL, method = "average")
```

## Arguments

- tcr_df:

  Data.frame with TCR columns.

- organism:

  Character string (`"human"` or `"mouse"`).

- k:

  Integer. Number of clusters. Exactly one of `k` or `h` must be
  specified.

- h:

  Numeric. Height at which to cut the dendrogram. Exactly one of `k` or
  `h` must be specified.

- method:

  Clustering method. Default `"average"`.

## Value

An integer vector of cluster assignments (length `nrow(tcr_df)`).

## See also

[`tcrdist_hclust`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_hclust.md),
[`neighborhood_test`](https://shihanli92.github.io/tcrdistR/reference/neighborhood_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
clusters <- cluster_tcrs(tcr_df, "human", k = 5)
table(clusters)
} # }
```
