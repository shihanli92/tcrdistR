# Sparse TCRdist matrix as COO triplets (C++ implementation)

Computes pairwise TCRdist distances for all (i, j) pairs with i \< j and
distance \<= `threshold`, returning a sparse representation in
coordinate (COO) format.

## Usage

``` r
rcpp_tcrdist_sparse(
  va_genes,
  cdr3a_seqs,
  vb_genes,
  cdr3b_seqs,
  v_dist_a,
  v_dist_b,
  threshold,
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

- threshold:

  Numeric. Maximum distance to include in the result. Pairs with
  distance \> threshold are omitted.

- weight_cdr3_region:

  Integer. CDR3 alignment distance multiplier. Default 3 matches
  `WEIGHT_CDR3_REGION`.

- gap_penalty_cdr3_region:

  Integer. Per-residue length-difference penalty. Default 12 matches
  `GAP_PENALTY_CDR3_REGION`.

## Value

A `List` with four elements:

- `i`:

  Integer vector of 1-based row indices.

- `j`:

  Integer vector of 1-based column indices.

- `x`:

  Numeric vector of distance values.

- `n`:

  Integer scalar: matrix dimension (N).

## Details

Three-stage early termination is used to skip pairs that cannot satisfy
the threshold:

1.  If V-region distance alone exceeds threshold, skip.

2.  If V-region + CDR3-alpha distance exceeds threshold, skip.

3.  If full distance exceeds threshold, skip.

`Rcpp::checkUserInterrupt()` is called every 100 outer-loop rows.

## Examples

``` r
if (FALSE) { # \dontrun{
  triplets <- rcpp_tcrdist_sparse(
    va_genes = c("TRAV1-1*01", "TRAV1-2*01"),
    cdr3a_seqs = c("CAVSANSGTYF", "CAVSANSGTYF"),
    vb_genes = c("TRBV20-1*01", "TRBV20-1*01"),
    cdr3b_seqs = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
    v_dist_a = v_alpha_mat, v_dist_b = v_beta_mat,
    threshold = 50
  )
} # }
```
