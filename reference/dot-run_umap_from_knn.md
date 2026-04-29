# Run UMAP from precomputed K-nearest-neighbor data

Passes precomputed neighbor indices and distances to
[`umap`](https://jlmelville.github.io/uwot/reference/umap.html) via its
`nn_method` list interface. Port of rconga's `.run_umap_from_knn()`,
adapted for 1-based indices.

## Usage

``` r
.run_umap_from_knn(
  knn_indices,
  knn_distances,
  n_components = 2L,
  min_dist = 0.5,
  spread = 1,
  seed = NULL,
  n_threads = 1L
)
```

## Arguments

- knn_indices:

  Integer matrix (N x K). 1-based neighbor indices.

- knn_distances:

  Numeric matrix (N x K). Distances.

- n_components:

  Integer. UMAP output dimensions. Default `2L`.

- min_dist:

  Numeric. UMAP min_dist. Default `0.5`.

- spread:

  Numeric. UMAP spread. Default `1.0`.

- seed:

  Integer or `NULL`. Random seed.

- n_threads:

  Integer. Threads for optimization. Default `1L`.

## Value

Numeric matrix (N x n_components).
