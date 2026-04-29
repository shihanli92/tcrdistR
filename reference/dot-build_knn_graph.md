# Build a fuzzy simplicial set graph from KNN data

Constructs a sparse N x N adjacency matrix using UMAP's fuzzy simplicial
set algorithm. Edge weights are `exp(-(d - rho) / sigma)` and the graph
is symmetrized via fuzzy union: `W + W^T - W * W^T`. Port of rconga's
`.build_knn_graph()`, adapted for 1-based indices.

## Usage

``` r
.build_knn_graph(knn_indices, knn_distances, n_cells)
```

## Arguments

- knn_indices:

  Integer matrix (N x K). 1-based neighbor indices (as returned by
  [`tcrdist_knn`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_knn.md)).

- knn_distances:

  Numeric matrix (N x K). Corresponding distances.

- n_cells:

  Integer scalar. Total number of observations.

## Value

Sparse `dgCMatrix` (N x N) with symmetrized fuzzy membership weights.
