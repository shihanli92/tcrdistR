# Plot a 2D scatter of TCR embeddings

Generic scatter plot for kernel PCA, UMAP, or any 2D embedding of TCR
repertoire data. Supports continuous (viridis) and categorical
(tab10/tab20) coloring with optional centroid labels, faceting, and
clone highlighting.

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
  na_color = "#DDDDDD",
  metadata = NULL,
  facet_by = NULL,
  highlight = NULL,
  highlight_color = "#E41A1C",
  background_alpha = 0.15
)
```

## Arguments

- coords:

  Numeric matrix with 2 columns (embedding coordinates).

- color_by:

  Optional vector of length `nrow(coords)`. Numeric for continuous
  coloring (viridis), factor/character for categorical coloring
  (tab10/tab20). `NULL` plots all points in gray. When `metadata` is
  provided, can be a column name string.

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

- metadata:

  Optional data.frame with `nrow(coords)` rows. When provided,
  `color_by` and `facet_by` can be column name strings, and `highlight`
  can be a named list of filter conditions applied to these columns.

- facet_by:

  Optional faceting variable(s). Either a character/factor vector for
  single-variable faceting (`facet_wrap`), or a data.frame with 1–2
  columns for grid faceting (`facet_grid`). When `metadata` is provided,
  can be a character vector of 1–2 column names.

- highlight:

  Optional. Which points to highlight, with non-highlighted points faded
  to a gray background. Accepts a logical vector, integer indices, or a
  named list for multi-column filtering (requires `metadata`; e.g.,
  `list(epitope = "PA", subject = c("S1", "S2"))`).

- highlight_color:

  Character. Color for highlighted points when `color_by` is `NULL`.
  Default `"#E41A1C"`.

- background_alpha:

  Numeric. Alpha for non-highlighted points. Default `0.15`.

## Value

A `ggplot` object.

## See also

[`compute_tcrdist_kernel_pca`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_kernel_pca.md),
[`compute_tcrdist_umap`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_umap.md),
[`plot_tcrdist_heatmap`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_heatmap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
kpca <- compute_tcrdist_kernel_pca(tcr_df, "human")
plot_tcr_scatter(kpca$embeddings[, 1:2], color_by = tcr_df$epitope)

# Faceting by epitope
plot_tcr_scatter(kpca$embeddings[, 1:2], color_by = tcr_df$epitope,
                 facet_by = tcr_df$epitope)

# Highlight specific clones
plot_tcr_scatter(kpca$embeddings[, 1:2],
                 highlight = tcr_df$epitope == "PA")

# Using metadata for column-name lookups
plot_tcr_scatter(kpca$embeddings[, 1:2],
                 metadata = tcr_df,
                 color_by = "epitope",
                 facet_by = c("epitope", "subject"),
                 highlight = list(epitope = "PA"))
} # }
```
