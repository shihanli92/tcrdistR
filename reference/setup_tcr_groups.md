# Assign alpha and beta chain group indices

For each TCR in `tcr_df`, builds a composite key from the chain gene
usage and CDR3 sequence, then assigns a 0-based group index so that
identical chains share the same group. These groups are used for
same-chain masking during neighbor search and Poisson testing.

## Usage

``` r
setup_tcr_groups(tcr_df)
```

## Arguments

- tcr_df:

  A `data.frame` with at least the columns `va`, `ja`, `cdr3a`, `vb`,
  `jb`, `cdr3b`. Optional columns `cdr3a_nucseq`, `cdr3b_nucseq`, and
  `subject_id` are included in the key when present.

## Value

A list with two elements:

- `agroups`:

  Integer vector of length `nrow(tcr_df)`. 0-based alpha-chain group
  indices.

- `bgroups`:

  Integer vector of length `nrow(tcr_df)`. 0-based beta-chain group
  indices.

## See also

[`find_clumping`](https://shihanli92.github.io/tcrdistR/reference/find_clumping.md)

## Examples

``` r
# \donttest{
tcr_df <- data.frame(
  va = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
  ja = c("TRAJ33*01", "TRAJ33*01", "TRAJ49*01"),
  cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
  vb = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
  jb = c("TRBJ2-7*01", "TRBJ2-7*01", "TRBJ1-1*01"),
  cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
  stringsAsFactors = FALSE
)
groups <- setup_tcr_groups(tcr_df)
groups$agroups  # c(0, 0, 1)  (first two share alpha chain)
#> [1] 0 0 1
# }
```
