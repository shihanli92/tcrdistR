# Build a TCR clone network from distance data

Constructs an igraph graph where vertices are TCR clones and edges
connect clones within a distance threshold. Edge weights encode
similarity via exponential decay: `exp(-distance / scale)`.

## Usage

``` r
compute_tcr_network(
  tcrs,
  organism = NULL,
  threshold = NULL,
  dist_matrix = NULL,
  scale = NULL,
  min_edges = 0L,
  jitter = TRUE,
  layout = "fr",
  seed = NULL
)
```

## Arguments

- tcrs:

  Data.frame with TCR columns (`va`, `cdr3a`, `vb`, `cdr3b`). All
  columns are attached as vertex attributes on the output graph.

- organism:

  Character string (`"human"` or `"mouse"`). Ignored when `dist_matrix`
  is provided.

- threshold:

  Numeric or `NULL`. Maximum TCRdist for an edge. If `NULL`,
  auto-detected from the distance distribution valley.

- dist_matrix:

  Optional precomputed N x N distance matrix (e.g. from
  [`tcrdist_matrix()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)
  or `TCRrep@paired_dist`). If provided, `organism` is ignored and
  distances are taken from this matrix.

- scale:

  Numeric or `NULL`. Scale parameter for the similarity transform
  `exp(-dist / scale)`. If `NULL`, defaults to `threshold / 4`.

- min_edges:

  Integer. Minimum number of edges a vertex must have to be kept. `0`
  (default) keeps all vertices. `1` removes singletons (isolated nodes),
  `2` removes vertices with fewer than 2 edges, etc.

- jitter:

  Logical. If `TRUE` (default), slightly offset overlapping nodes (e.g.
  identical clones with distance 0) so they are individually visible
  instead of stacking on top of each other.

- layout:

  Character. Layout algorithm: `"fr"` (Fruchterman-Reingold), `"kk"`
  (Kamada-Kawai), `"drl"`, `"circle"`, or `"grid"`. Default `"fr"`.

- seed:

  Integer or `NULL`. Random seed for layout reproducibility.

## Value

A named list with elements:

- `graph`:

  An `igraph` object. Vertex attributes include all columns from `tcrs`.
  Edge attributes: `weight` (similarity) and `distance` (TCRdist).

- `layout`:

  Numeric matrix (N x 2) of layout coordinates.

- `threshold`:

  Numeric. Threshold used.

- `scale`:

  Numeric. Scale parameter used.

- `n_components`:

  Integer. Number of connected components.

- `n_edges`:

  Integer. Number of edges in the graph.

- `dist_plot`:

  `ggplot` object showing the distance distribution with the threshold
  line, or `NULL` if ggplot2 is not available.

## Details

When `threshold = NULL`, the threshold is auto-detected by finding the
valley between the two peaks of the (typically bimodal) TCRdist
distribution, using kernel density estimation on a subsample.

A distance distribution plot with the threshold marked is automatically
displayed when ggplot2 is available. The plot is also stored in the
returned list as `dist_plot`.

## See also

[`plot_tcr_network`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_network.md),
[`tcrdist_sparse`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_sparse.md),
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)

## Examples

``` r
# \donttest{
data(dash)
sub <- dash[1:200, ]
net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42)

net$n_components
#> [1] 146
net$n_edges
#> [1] 108
# }
```
