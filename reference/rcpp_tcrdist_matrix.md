# Pairwise TCRdist matrix (C++ implementation)

Computes the full N x N symmetric matrix of paired-chain TCRdist
distances. For each pair `(i, j)`: \$\$ d(i,j) = v\\dist\\a\[va_i,
va_j\] + cdr3\\dist(cdr3a_i, cdr3a_j) + v\\dist\\b\[vb_i, vb_j\] +
cdr3\\dist(cdr3b_i, cdr3b_j) \$\$

## Usage

``` r
rcpp_tcrdist_matrix(
  va_genes,
  cdr3a_seqs,
  vb_genes,
  cdr3b_seqs,
  v_dist_a,
  v_dist_b,
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
  for alpha genes (row/column names are V-gene allele strings, e.g.
  `"TRAV1-1*01"`).

- v_dist_b:

  Named square `NumericMatrix`. Same structure for beta genes.

- weight_cdr3_region:

  Integer. CDR3 alignment distance multiplier. Default 3 matches
  `WEIGHT_CDR3_REGION`.

- gap_penalty_cdr3_region:

  Integer. Per-residue length-difference penalty. Default 12 matches
  `GAP_PENALTY_CDR3_REGION`.

## Value

An N x N symmetric `NumericMatrix` of TCRdist distances.

## Details

Only the upper triangle is computed; it is mirrored to the lower
triangle. The diagonal is zero. `Rcpp::checkUserInterrupt()` is called
every 100 outer-loop rows to allow the user to interrupt long
computations.

BSD4 values are compiled in as `constexpr`; no `bsd4` matrix argument is
needed (unlike the rconga equivalent).

## Examples

``` r
if (FALSE) { # \dontrun{
  mat <- rcpp_tcrdist_matrix(
    va_genes = c("TRAV1-1*01", "TRAV1-2*01"),
    cdr3a_seqs = c("CAVSANSGTYF", "CAVSANSGTYF"),
    vb_genes = c("TRBV20-1*01", "TRBV20-1*01"),
    cdr3b_seqs = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
    v_dist_a = v_alpha_mat,
    v_dist_b = v_beta_mat
  )
} # }
```
