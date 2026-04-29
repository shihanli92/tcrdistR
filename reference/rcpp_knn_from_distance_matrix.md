# K-nearest-neighbors from a precomputed distance matrix (C++ implementation)

For each row of a square distance matrix `D`, finds the K nearest
neighbors while masking out entries that share the same alpha-chain or
beta-chain group assignment. Same-group pairs are assigned the sentinel
distance 1e3 before selection.

## Usage

``` r
rcpp_knn_from_distance_matrix(D, K, agroups, bgroups, sort_nbrs = TRUE)
```

## Arguments

- D:

  Numeric matrix (N x N). Precomputed pairwise distance matrix. Must be
  square.

- K:

  Integer. Number of nearest neighbors to extract per row.

- agroups:

  Integer vector of length N. Alpha-chain group assignments. Rows
  sharing the same agroups value are masked from each other.

- bgroups:

  Integer vector of length N. Beta-chain group assignments. Rows sharing
  the same bgroups value are masked from each other.

- sort_nbrs:

  Logical. If `TRUE`, sort the K neighbors by ascending distance.
  Default `TRUE`.

## Value

A `List` with two elements:

- `knn_indices`:

  Integer matrix (N x K). 1-based neighbor indices.

- `knn_distances`:

  Numeric matrix (N x K). Corresponding distances.

## Details

Uses `std::nth_element` for O(N) partial sorting per row, optionally
followed by `std::sort` of the K selected neighbors.

## Examples

``` r
if (FALSE) { # \dontrun{
  D <- matrix(c(0,1,2,1,0,3,2,3,0), nrow=3)
  result <- rcpp_knn_from_distance_matrix(D, K=1L, agroups=1:3, bgroups=1:3)
} # }
```
