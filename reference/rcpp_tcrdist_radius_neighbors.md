# Radius-based TCRdist neighbor search with group masking (C++ implementation)

For each of the N input TCRs, finds all other TCRs within `radius`
TCRdist distance, excluding TCRs that share the same alpha or beta group
(same-group exclusion) and self-pairs.

## Usage

``` r
rcpp_tcrdist_radius_neighbors(
  va_genes,
  cdr3a_seqs,
  vb_genes,
  cdr3b_seqs,
  v_dist_a,
  v_dist_b,
  radius,
  agroups,
  bgroups,
  weight_cdr3_a = 3L,
  gap_penalty_cdr3_a = 12L,
  weight_cdr3_b = 3L,
  gap_penalty_cdr3_b = 12L
)
```

## Arguments

- va_genes:

  Character vector of length N. Alpha V-gene allele names.

- cdr3a_seqs:

  Character vector of length N. Alpha CDR3 sequences.

- vb_genes:

  Character vector of length N. Beta V-gene allele names.

- cdr3b_seqs:

  Character vector of length N. Beta CDR3 sequences.

- v_dist_a:

  Named square `NumericMatrix`. Pre-computed pairwise V-region distances
  for alpha genes.

- v_dist_b:

  Named square `NumericMatrix`. Pre-computed pairwise V-region distances
  for beta genes.

- radius:

  Numeric. Search radius. Only neighbors within this distance
  (inclusive) are returned.

- agroups:

  Integer vector of length N. Alpha-chain group assignments.

- bgroups:

  Integer vector of length N. Beta-chain group assignments.

- weight_cdr3_region:

  Integer. CDR3 alignment distance multiplier. Default 3 matches
  `WEIGHT_CDR3_REGION`.

- gap_penalty_cdr3_region:

  Integer. Per-residue length-difference penalty. Default 12 matches
  `GAP_PENALTY_CDR3_REGION`.

## Value

A `List` of length N. Each element is a `List` with:

- `indices`:

  Integer vector of 1-based neighbor indices.

- `distances`:

  Numeric vector of corresponding distances.

## Details

Three-stage early termination is applied per pair:

1.  If V-region distance alone exceeds `radius`, skip.

2.  If V-region + CDR3-alpha distance exceeds `radius`, skip.

3.  If full distance exceeds `radius`, skip.

`Rcpp::checkUserInterrupt()` is called every 100 rows.

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- rcpp_tcrdist_radius_neighbors(
    va_genes = c("TRAV1-1*01", "TRAV1-2*01"),
    cdr3a_seqs = c("CAVSANSGTYF", "CAVSANSGTYF"),
    vb_genes = c("TRBV20-1*01", "TRBV20-1*01"),
    cdr3b_seqs = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
    v_dist_a = v_alpha_mat, v_dist_b = v_beta_mat,
    radius = 50, agroups = 1:2, bgroups = 1:2
  )
} # }
```
