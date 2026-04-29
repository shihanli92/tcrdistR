# K-nearest-neighbors from a precomputed distance matrix

Extracts K nearest neighbors for each row of a precomputed N x N
distance matrix `D`, with optional same-group masking.

## Usage

``` r
knn_from_matrix(D, K, agroups = NULL, bgroups = NULL, sort_nbrs = TRUE)
```

## Arguments

- D:

  Numeric matrix (N x N). Precomputed pairwise distance matrix. Must be
  square.

- K:

  Integer. Number of nearest neighbors to return. Must satisfy
  `1 <= K <= N - 1`.

- agroups:

  Integer vector of length N, or `NULL`. Alpha-chain group assignments.
  `NULL` assigns unique groups (no masking).

- bgroups:

  Integer vector of length N, or `NULL`. Beta-chain group assignments.

- sort_nbrs:

  Logical. If `TRUE` (default), sort K neighbors by ascending distance.

## Value

A `list` with two elements:

- `knn_indices`:

  Integer matrix (N x K). 1-based indices.

- `knn_distances`:

  Numeric matrix (N x K). Distances.

## See also

[`knn_from_pca`](https://shihanli92.github.io/tcrdistR/reference/knn_from_pca.md),
[`tcrdist_knn`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_knn.md)

## Examples

``` r
# \donttest{
tcrs <- data.frame(
  va    = c("TRAV1-1*01", "TRAV1-1*01"),
  cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
  vb    = c("TRBV19*01", "TRBV19*01"),
  cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
  stringsAsFactors = FALSE
)
D   <- tcrdist_matrix(tcrs, "human")
knn <- knn_from_matrix(D, K = 1L)
# }
```
