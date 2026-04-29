# K-nearest-neighbors from a PCA embedding matrix

Computes K nearest neighbors by Euclidean distance from an N x D PCA (or
other embedding) matrix, without constructing the full N x N distance
matrix. Optional same-group masking applies.

## Usage

``` r
knn_from_pca(pca_matrix, K, agroups = NULL, bgroups = NULL, sort_nbrs = TRUE)
```

## Arguments

- pca_matrix:

  Numeric matrix (N x D). Rows are samples, columns are embedding
  dimensions.

- K:

  Integer. Number of nearest neighbors. Must satisfy `1 <= K <= N - 1`.

- agroups:

  Integer vector of length N, or `NULL`. `NULL` assigns unique groups
  (no masking).

- bgroups:

  Integer vector of length N, or `NULL`.

- sort_nbrs:

  Logical. If `TRUE` (default), sort K neighbors by ascending Euclidean
  distance.

## Value

A `list` with two elements:

- `knn_indices`:

  Integer matrix (N x K). 1-based indices.

- `knn_distances`:

  Numeric matrix (N x K). Euclidean distances.

## See also

[`knn_from_matrix`](https://shihanli92.github.io/tcrdistR/reference/knn_from_matrix.md),
[`compute_tcrdist_kernel_pca`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_kernel_pca.md)

## Examples

``` r
# \donttest{
pca <- matrix(rnorm(30), nrow = 10, ncol = 3)
knn <- knn_from_pca(pca, K = 3L)
# }
```
