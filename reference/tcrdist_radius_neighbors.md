# Radius-based TCRdist neighbor search with optional group masking

For each of the N input TCRs, finds all other TCRs within `radius`
TCRdist distance. TCRs sharing the same `agroups` or `bgroups` value are
excluded (same-group masking).

## Usage

``` r
tcrdist_radius_neighbors(
  tcrs,
  organism,
  radius,
  agroups = NULL,
  bgroups = NULL,
  components = "all",
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

- radius:

  Numeric. Search radius (inclusive). Pairs with distance \<= `radius`
  are returned. Must be \>= 0.

- agroups:

  Integer vector of length N, or `NULL`. Alpha-chain group assignments.
  `NULL` assigns each TCR a unique group (no masking).

- bgroups:

  Integer vector of length N, or `NULL`. Beta-chain group assignments.

- components:

  Character. Which distance components to include. See
  [`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)
  for details. Default `"all"`.

- weight_cdr3:

  Integer. CDR3 distance weight. Defaults to `WEIGHT_CDR3_REGION` (3L).

- gap_penalty_cdr3:

  Integer. CDR3 gap penalty. Defaults to `GAP_PENALTY_CDR3_REGION`
  (12L).

## Value

A `list` of length N. Each element `[[i]]` is a `list` with:

- `indices`:

  Integer vector. 1-based indices of neighbors within `radius`.

- `distances`:

  Numeric vector. Corresponding distances.

## See also

[`tcrdist_knn`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_knn.md),
[`find_clumping`](https://shihanli92.github.io/tcrdistR/reference/find_clumping.md)

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
nbrs <- tcrdist_radius_neighbors(tcrs, "human", radius = 50)
# }
```
