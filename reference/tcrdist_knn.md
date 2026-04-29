# K-nearest-neighbors by TCRdist with optional group masking

For each of the N input TCRs, finds the K nearest neighbors by TCRdist
distance. TCRs sharing the same `agroups` or `bgroups` value are
excluded from each other's neighborhood (same-group masking). When
`agroups` and `bgroups` are `NULL` (default), every TCR has its own
unique group so no masking occurs.

## Usage

``` r
tcrdist_knn(
  tcrs,
  organism,
  K,
  agroups = NULL,
  bgroups = NULL,
  sort_nbrs = TRUE,
  weight_cdr3 = WEIGHT_CDR3_REGION,
  gap_penalty_cdr3 = GAP_PENALTY_CDR3_REGION
)
```

## Arguments

- tcrs:

  A `data.frame` with at least the following columns:

  `va`

  :   Character. Alpha-chain V-gene allele.

  `cdr3a`

  :   Character. Alpha-chain CDR3 amino acid sequence.

  `vb`

  :   Character. Beta-chain V-gene allele.

  `cdr3b`

  :   Character. Beta-chain CDR3 amino acid sequence.

- organism:

  Character string. Organism key, e.g. `"human"`.

- K:

  Integer. Number of nearest neighbors to return per TCR. Must satisfy
  `1 <= K <= N - 1`.

- agroups:

  Integer vector of length N, or `NULL`. Alpha-chain group assignments.
  TCRs with the same value are masked from each other. `NULL` assigns
  each TCR a unique group (no masking).

- bgroups:

  Integer vector of length N, or `NULL`. Beta-chain group assignments.
  Same semantics as `agroups`.

- sort_nbrs:

  Logical. If `TRUE` (default), sort each row's K neighbors by ascending
  distance.

- weight_cdr3:

  Integer. CDR3 distance weight. Defaults to `WEIGHT_CDR3_REGION` (3L).

- gap_penalty_cdr3:

  Integer. CDR3 gap penalty. Defaults to `GAP_PENALTY_CDR3_REGION`
  (12L).

## Value

A `list` with two elements:

- `knn_indices`:

  Integer matrix (N x K). 1-based row indices of the K nearest neighbors
  for each TCR.

- `knn_distances`:

  Numeric matrix (N x K). Corresponding TCRdist distances.

## See also

[`tcrdist_radius_neighbors`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_radius_neighbors.md),
[`knn_from_matrix`](https://shihanli92.github.io/tcrdistR/reference/knn_from_matrix.md),
[`knn_from_pca`](https://shihanli92.github.io/tcrdistR/reference/knn_from_pca.md)

## Examples

``` r
# \donttest{
tcrs <- data.frame(
  va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
  cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
  vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
  cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
  stringsAsFactors = FALSE
)
knn <- tcrdist_knn(tcrs, "human", K = 1L)
# }
```
