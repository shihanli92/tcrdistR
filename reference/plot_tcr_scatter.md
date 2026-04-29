# Plot a 2D scatter of TCR embeddings

Generic scatter plot for kernel PCA, UMAP, or any 2D embedding of TCR
repertoire data. Supports continuous (viridis) and categorical
(tab10/tab20) coloring with optional centroid labels.

## Usage

``` r
plot_tcr_scatter(
  coords,
  color_by = NULL,
  title = NULL,
  point_size = 1,
  alpha = 1,
  axis_label_prefix = "KPCA",
  legend_title = NULL,
  palette = NULL,
  show_labels = FALSE,
  label_size = 3,
  na_color = "#DDDDDD"
)
```

## Arguments

- coords:

  Numeric matrix with 2 columns (embedding coordinates).

- color_by:

  Optional vector of length `nrow(coords)`. Numeric for continuous
  coloring (viridis), factor/character for categorical coloring
  (tab10/tab20). `NULL` plots all points in gray.

- title:

  Optional plot title.

- point_size:

  Numeric. Point size. Default `1`.

- alpha:

  Numeric. Point opacity. Default `1`.

- axis_label_prefix:

  Character string. Prefix for axis labels. Default `"KPCA"`.

- legend_title:

  Optional legend title.

- palette:

  Character vector of colors to override default palette.

- show_labels:

  Logical. If `TRUE` and `color_by` is categorical, add centroid labels.
  Default `FALSE`.

- label_size:

  Numeric. Label text size. Default `3`.

- na_color:

  Character. Color for NA values. Default `"#DDDDDD"`.

## Value

A `ggplot` object.

## See also

[`compute_tcrdist_kernel_pca`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_kernel_pca.md),
[`plot_tcrdist_heatmap`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_heatmap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
kpca <- compute_tcrdist_kernel_pca(tcr_df, "human")
plot_tcr_scatter(kpca$embeddings[, 1:2], color_by = tcr_df$epitope)
} # }
```
