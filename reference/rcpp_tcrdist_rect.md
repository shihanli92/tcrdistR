# Rectangular TCRdist query-vs-reference distance matrix (C++ implementation)

Computes an nq x nr matrix of paired-chain TCRdist distances between a
query set of `nq` TCRs and a reference set of `nr` TCRs. For each
query-reference pair `(i, j)`: \$\$ d(i,j) = v\\dist\\a\[va\\query_i,
va\\ref_j\] + cdr3\\dist(cdr3a\\query_i, cdr3a\\ref_j) +
v\\dist\\b\[vb\\query_i, vb\\ref_j\] + cdr3\\dist(cdr3b\\query_i,
cdr3b\\ref_j) \$\$

## Usage

``` r
rcpp_tcrdist_rect(
  query_va,
  query_cdr3a,
  query_vb,
  query_cdr3b,
  ref_va,
  ref_cdr3a,
  ref_vb,
  ref_cdr3b,
  v_dist_a,
  v_dist_b,
  weight_cdr3_a = 3L,
  gap_penalty_cdr3_a = 12L,
  weight_cdr3_b = 3L,
  gap_penalty_cdr3_b = 12L
)
```

## Arguments

- query_va:

  Character vector of length nq. Query alpha V-gene alleles.

- query_cdr3a:

  Character vector of length nq. Query alpha CDR3 sequences.

- query_vb:

  Character vector of length nq. Query beta V-gene alleles.

- query_cdr3b:

  Character vector of length nq. Query beta CDR3 sequences.

- ref_va:

  Character vector of length nr. Reference alpha V-gene alleles.

- ref_cdr3a:

  Character vector of length nr. Reference alpha CDR3 sequences.

- ref_vb:

  Character vector of length nr. Reference beta V-gene alleles.

- ref_cdr3b:

  Character vector of length nr. Reference beta CDR3 sequences.

- v_dist_a:

  Named square `NumericMatrix`. Pre-computed pairwise V-region distances
  for alpha genes.

- v_dist_b:

  Named square `NumericMatrix`. Pre-computed pairwise V-region distances
  for beta genes.

- weight_cdr3_region:

  Integer. CDR3 alignment distance multiplier. Default 3 matches
  `WEIGHT_CDR3_REGION`.

- gap_penalty_cdr3_region:

  Integer. Per-residue length-difference penalty. Default 12 matches
  `GAP_PENALTY_CDR3_REGION`.

## Value

An nq x nr `NumericMatrix` of TCRdist distances.

## Details

No symmetry is exploited; all nq \* nr distances are computed.
`Rcpp::checkUserInterrupt()` is called every 100 outer-loop rows.

## Examples

``` r
if (FALSE) { # \dontrun{
  mat <- rcpp_tcrdist_rect(
    query_va = c("TRAV1-1*01"),
    query_cdr3a = c("CAVSANSGTYF"),
    query_vb = c("TRBV20-1*01"),
    query_cdr3b = c("CASSIRSSYEQYF"),
    ref_va = c("TRAV1-1*01", "TRAV1-2*01"),
    ref_cdr3a = c("CAVSANSGTYF", "CAVSANSGTYF"),
    ref_vb = c("TRBV20-1*01", "TRBV20-1*01"),
    ref_cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
    v_dist_a = v_alpha_mat,
    v_dist_b = v_beta_mat
  )
} # }
```
