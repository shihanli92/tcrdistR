# Plot the distribution of pairwise TCRdist distances

Histogram with overlaid density curve of the upper triangle of a
pairwise distance matrix. A vertical dashed line marks the median
distance. When `threshold` is provided, a red vertical line is drawn at
that value with a label.

## Usage

``` r
plot_distance_distribution(
  dist_matrix,
  title = NULL,
  binwidth = NULL,
  threshold = NULL
)
```

## Arguments

- dist_matrix:

  Numeric matrix or numeric vector. If a square matrix, the upper
  triangle is extracted; if a vector, used directly.

- title:

  Optional plot title.

- binwidth:

  Numeric. Histogram bin width. If `NULL`, uses ggplot2 default.

- threshold:

  Numeric or `NULL`. If not `NULL`, a vertical line is drawn at this
  distance value and labeled.

## Value

A `ggplot` object.

## See also

[`plot_tcrdist_heatmap`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_heatmap.md),
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md),
[`compute_tcr_network`](https://shihanli92.github.io/tcrdistR/reference/compute_tcr_network.md)

## Examples

``` r
if (FALSE) { # \dontrun{
mat <- matrix(c(0,10,20, 10,0,15, 20,15,0), nrow = 3)
plot_distance_distribution(mat, threshold = 12)
} # }
```
