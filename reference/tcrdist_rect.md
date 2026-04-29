# Compute rectangular TCRdist matrix between a query set and a reference set

Computes an nq x nr matrix of paired-chain TCRdist distances between a
query set of `nq` TCRs and a reference set of `nr` TCRs. Unlike
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md),
the result is not symmetric because query and reference need not be the
same set.

## Usage

``` r
tcrdist_rect(
  query,
  ref,
  organism,
  weight_cdr3 = WEIGHT_CDR3_REGION,
  gap_penalty_cdr3 = GAP_PENALTY_CDR3_REGION
)
```

## Arguments

- query:

  A `data.frame` with at least the following columns:

  `va`

  :   Character. Alpha-chain V-gene allele, e.g. `"TRAV1-1*01"`.

  `cdr3a`

  :   Character. Alpha-chain CDR3 amino acid sequence.

  `vb`

  :   Character. Beta-chain V-gene allele, e.g. `"TRBV19*01"`.

  `cdr3b`

  :   Character. Beta-chain CDR3 amino acid sequence.

- ref:

  A `data.frame` with the same four columns as `query`.

- organism:

  Character string. Organism key understood by `load_gene_database`,
  e.g. `"human"` or `"mouse"`.

- weight_cdr3:

  Integer. Weight applied to CDR3 distances. Defaults to
  `WEIGHT_CDR3_REGION` (3L).

- gap_penalty_cdr3:

  Integer. Gap penalty for CDR3 alignments. Defaults to
  `GAP_PENALTY_CDR3_REGION` (12L).

## Value

A numeric matrix of dimensions nq x nr. Row names are row indices of
`query` as character strings; column names are row indices of `ref`.

## Details

For each query-reference pair `(i, j)` the distance is: \$\$ d(i,j) =
v\\dist\\a\[\mathrm{va\\query}\_i, \mathrm{va\\ref}\_j\] +
\mathrm{cdr3\\dist}(\mathrm{cdr3a\\query}\_i, \mathrm{cdr3a\\ref}\_j) +
v\\dist\\b\[\mathrm{vb\\query}\_i, \mathrm{vb\\ref}\_j\] +
\mathrm{cdr3\\dist}(\mathrm{cdr3b\\query}\_i, \mathrm{cdr3b\\ref}\_j)
\$\$

When `query` and `ref` are the same data.frame the result equals
`tcrdist_matrix(query, organism)` (no symmetry short-cut is taken so the
diagonal may differ by floating-point rounding from zero, but is
numerically identical when both data.frames are identical objects).

## See also

[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md),
[`tcrdist_sparse`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_sparse.md),
[`tcrdist_join`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_join.md)

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
# Self-comparison == tcrdist_matrix
rect <- tcrdist_rect(tcrs, tcrs, "human")
mat  <- tcrdist_matrix(tcrs, "human")
# }
```
