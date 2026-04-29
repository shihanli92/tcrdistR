# K-nearest-neighbors by TCRdist with group masking (C++ implementation)

For each of the N input TCRs, finds the K nearest neighbors by TCRdist
while masking out TCRs that share the same alpha-group or beta-group
assignment (same-group exclusion). Same-group pairs and self-pairs are
assigned the sentinel distance 1e3 before selection.

## Usage

``` r
rcpp_tcrdist_knn(
  va_genes,
  cdr3a_seqs,
  vb_genes,
  cdr3b_seqs,
  v_dist_a,
  v_dist_b,
  K,
  agroups,
  bgroups,
  sort_nbrs = TRUE,
  weight_cdr3_region = 3L,
  gap_penalty_cdr3_region = 12L
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

- K:

  Integer. Number of nearest neighbors to return per TCR.

- agroups:

  Integer vector of length N. Alpha-chain group assignments. TCRs
  sharing the same agroups value are masked from each other.

- bgroups:

  Integer vector of length N. Beta-chain group assignments. TCRs sharing
  the same bgroups value are masked from each other.

- sort_nbrs:

  Logical. If `TRUE`, sort the K neighbors by ascending distance.
  Default `TRUE`.

- weight_cdr3_region:

  Integer. CDR3 alignment distance multiplier. Default 3 matches
  `WEIGHT_CDR3_REGION`.

- gap_penalty_cdr3_region:

  Integer. Per-residue length-difference penalty. Default 12 matches
  `GAP_PENALTY_CDR3_REGION`.

## Value

A `List` with two elements:

- `knn_indices`:

  Integer matrix (N x K). 1-based neighbor indices.

- `knn_distances`:

  Numeric matrix (N x K). Corresponding distances.

## Details

Uses `std::nth_element` for O(N) partial sorting per row, optionally
followed by O(K log K) full sorting of the K selected neighbors.
`Rcpp::checkUserInterrupt()` is called every 100 rows.

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- rcpp_tcrdist_knn(
    va_genes = c("TRAV1-1*01", "TRAV1-2*01"),
    cdr3a_seqs = c("CAVSANSGTYF", "CAVSANSGTYF"),
    vb_genes = c("TRBV20-1*01", "TRBV20-1*01"),
    cdr3b_seqs = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
    v_dist_a = v_alpha_mat, v_dist_b = v_beta_mat,
    K = 1L,
    agroups = 1:2, bgroups = 1:2
  )
} # }
```
