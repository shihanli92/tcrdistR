# Plot a TCRdist pairwise distance heatmap

Renders a heatmap of a pairwise distance matrix using
[`ggplot2::geom_tile()`](https://ggplot2.tidyverse.org/reference/geom_tile.html).
Optionally clusters rows and columns via hierarchical clustering
(UPGMA).

## Usage

``` r
plot_tcrdist_heatmap(dist_matrix, labels = NULL, cluster = TRUE, title = NULL)
```

## Arguments

- dist_matrix:

  Numeric matrix. Square pairwise distance matrix.

- labels:

  Character vector. Row/column labels. If `NULL`, uses
  `rownames(dist_matrix)` or integer indices.

- cluster:

  Logical. If `TRUE` (default), reorder rows and columns by hierarchical
  clustering.

- title:

  Optional plot title.

## Value

A `ggplot` object.

## See also

[`plot_tcrdist_dendrogram`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_dendrogram.md),
[`plot_distance_distribution`](https://shihanli92.github.io/tcrdistR/reference/plot_distance_distribution.md),
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
mat <- matrix(c(0,10,20, 10,0,15, 20,15,0), nrow = 3)
plot_tcrdist_heatmap(mat, labels = c("TCR1", "TCR2", "TCR3"))
} # }
```
