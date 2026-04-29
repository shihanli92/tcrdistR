# Compute sparse TCRdist matrix as a `dgCMatrix`

Computes pairwise TCRdist distances for all pairs `(i, j)` where `i < j`
and the distance is at most `threshold`, then returns a symmetric sparse
matrix in compressed-column format (`dgCMatrix`).

## Usage

``` r
tcrdist_sparse(
  tcrs,
  organism,
  threshold,
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

- threshold:

  Numeric. Maximum distance to include. Pairs with distance strictly
  greater than `threshold` are omitted. Must be \>= 0.

- weight_cdr3:

  Integer. CDR3 distance weight. Defaults to `WEIGHT_CDR3_REGION` (3L).

- gap_penalty_cdr3:

  Integer. CDR3 gap penalty. Defaults to `GAP_PENALTY_CDR3_REGION`
  (12L).

## Value

A symmetric sparse matrix of class `dgCMatrix` with dimensions N x N.
Off-diagonal entries (i, j) and (j, i) are present for all pairs within
the threshold. The diagonal is structural zero (not stored). Returns an
all-zero N x N `dgCMatrix` if no pairs satisfy the threshold.

## Details

Three-stage early termination is used internally:

1.  If the V-region distance sum alone exceeds `threshold`, skip.

2.  If V-region + CDR3-alpha distance exceeds `threshold`, skip.

3.  If the full distance exceeds `threshold`, skip.

For `threshold = Inf` all pairs are evaluated (equivalent to a full
dense matrix stored as sparse), which serves as a correctness check
against
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md).

## See also

[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md),
[`tcrdist_rect`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_rect.md),
[`tcrdist_radius_neighbors`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_radius_neighbors.md)

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
sp <- tcrdist_sparse(tcrs, "human", threshold = 50)
# }
```
