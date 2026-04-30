# Compute pairwise TCRdist distance matrix

Computes the full N x N symmetric matrix of paired-chain TCRdist
distances for a collection of TCRs. Each off-diagonal entry `[i, j]`
equals the sum of V-alpha, CDR3-alpha, V-beta, and CDR3-beta component
distances between TCR `i` and TCR `j`. Computation is dispatched to a
C++ implementation for performance.

## Usage

``` r
tcrdist_matrix(
  tcrs,
  organism,
  components = "all",
  weight_cdr3 = WEIGHT_CDR3_REGION,
  gap_penalty_cdr3 = GAP_PENALTY_CDR3_REGION
)
```

## Arguments

- tcrs:

  A `data.frame` with at least the following columns:

  `va`

  :   Character. Alpha-chain V-gene allele, e.g. `"TRAV1-1*01"`.

  `cdr3a`

  :   Character. Alpha-chain CDR3 amino acid sequence.

  `vb`

  :   Character. Beta-chain V-gene allele, e.g. `"TRBV19*01"`.

  `cdr3b`

  :   Character. Beta-chain CDR3 amino acid sequence.

- organism:

  Character string. Organism key understood by `load_gene_database`,
  e.g. `"human"` or `"mouse"`.

- components:

  Character. Which distance components to include. Presets: `"all"`
  (default), `"cdr3"` (CDR3 only), `"v_region"` (CDR1+CDR2+CDR2.5 only),
  `"alpha"` (alpha chain), `"beta"` (beta chain). Or a character vector
  of individual terms: `"va"`, `"cdr3a"`, `"vb"`, `"cdr3b"`.

- weight_cdr3:

  Integer. Weight applied to CDR3 distances. Defaults to
  `WEIGHT_CDR3_REGION` (3L).

- gap_penalty_cdr3:

  Integer. Gap penalty for CDR3 alignments. Defaults to
  `GAP_PENALTY_CDR3_REGION` (12L).

## Value

A symmetric numeric matrix of dimensions N x N where N is `nrow(tcrs)`.
Row and column names are the row indices of `tcrs` as character strings.

## Details

The diagonal is zero (distance of a TCR to itself). The matrix is
symmetric by construction.

## See also

[`tcrdist_sparse`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_sparse.md),
[`tcrdist_rect`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_rect.md),
[`tcrdist_knn`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_knn.md),
[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md)

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
mat <- tcrdist_matrix(tcrs, "human")

# CDR3-only distance
mat_cdr3 <- tcrdist_matrix(tcrs, "human", components = "cdr3")

# Alpha chain only
mat_alpha <- tcrdist_matrix(tcrs, "human", components = "alpha")
# }
```
