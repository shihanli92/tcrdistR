# Plot the distribution of pairwise TCRdist distances

Histogram with overlaid density curve of the upper triangle of a
pairwise distance matrix. A vertical dashed line marks the median
distance.

## Usage

``` r
plot_distance_distribution(dist_matrix, title = NULL, binwidth = NULL)
```

## Arguments

- dist_matrix:

  Numeric matrix. Square pairwise distance matrix.

- title:

  Optional plot title.

- binwidth:

  Numeric. Histogram bin width. If `NULL`, uses ggplot2 default.

## Value

A `ggplot` object.

## See also

[`plot_tcrdist_heatmap`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_heatmap.md),
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
mat <- matrix(c(0,10,20, 10,0,15, 20,15,0), nrow = 3)
plot_distance_distribution(mat)
} # }
```
