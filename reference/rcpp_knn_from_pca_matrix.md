# K-nearest-neighbors from a PCA embedding matrix (C++ implementation)

Computes K nearest neighbors directly from PCA embeddings (N x D matrix)
without materializing the full N x N Euclidean distance matrix. For each
row, computes Euclidean distances to all N points on the fly, applies
group masking, and extracts the K nearest via `std::nth_element`.

## Usage

``` r
rcpp_knn_from_pca_matrix(pca_matrix, K, agroups, bgroups, sort_nbrs = TRUE)
```

## Arguments

- pca_matrix:

  Numeric matrix (N x D). PCA or other embedding coordinates. Rows are
  samples, columns are dimensions.

- K:

  Integer. Number of nearest neighbors to extract per point.

- agroups:

  Integer vector of length N. Alpha-chain group assignments.

- bgroups:

  Integer vector of length N. Beta-chain group assignments.

- sort_nbrs:

  Logical. If `TRUE`, sort the K neighbors by ascending Euclidean
  distance. Default `TRUE`.

## Value

A `List` with two elements:

- `knn_indices`:

  Integer matrix (N x K). 1-based neighbor indices.

- `knn_distances`:

  Numeric matrix (N x K). Euclidean distances.

## Details

Memory: O(N*D + N*K) instead of O(N^2).

## Examples

``` r
if (FALSE) { # \dontrun{
  pca <- matrix(rnorm(30), nrow=10, ncol=3)
  result <- rcpp_knn_from_pca_matrix(pca, K=3L, agroups=1:10, bgroups=1:10)
} # }
```
