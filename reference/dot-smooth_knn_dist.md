# Compute per-point bandwidth for fuzzy simplicial set

Binary search for sigma (per-point bandwidth) and rho (local
connectivity distance) used by UMAP's fuzzy simplicial set construction.
Port of rconga's `.smooth_knn_dist()`.

## Usage

``` r
.smooth_knn_dist(
  knn_distances,
  local_connectivity = 1,
  bandwidth = 1,
  n_iter = 64L
)
```

## Arguments

- knn_distances:

  Numeric matrix (N x K). Sorted KNN distances.

- local_connectivity:

  Numeric scalar. Default `1.0`.

- bandwidth:

  Numeric scalar. Default `1.0`.

- n_iter:

  Integer. Max binary search iterations. Default `64L`.

## Value

Named list with `rho` and `sigma` (numeric vectors of length N).
