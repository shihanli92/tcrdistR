# Find significant TCRdist matches between query and database TCRs

Computes paired TCRdist distances between query and database TCRs and
converts them to p-values adjusted for both the number of query and
database TCRs. Background distributions are estimated via the V(D)J
rearrangement model in
[`.estimate_background_tcrdist_distributions()`](https://shihanli92.github.io/tcrdistR/reference/dot-estimate_background_tcrdist_distributions.md).

## Usage

``` r
find_significant_tcrdist_matches(
  query_tcrs_df,
  db_tcrs_df,
  organism,
  adjusted_pvalue_threshold = 1,
  background_tcrs_df = NULL,
  num_random_samples = 50000L,
  fixup_alleles = TRUE
)
```

## Arguments

- query_tcrs_df:

  A data.frame with at least columns `va` (or `va_gene`), `cdr3a`, `vb`
  (or `vb_gene`), `cdr3b`. If `background_tcrs_df` is NULL, also needs
  `ja`, `jb`, `cdr3a_nucseq`, `cdr3b_nucseq`.

- db_tcrs_df:

  A data.frame with at least columns `va` (or `va_gene`), `cdr3a`, `vb`
  (or `vb_gene`), `cdr3b`.

- organism:

  Character string (e.g. `"human"`, `"mouse"`).

- adjusted_pvalue_threshold:

  Numeric. Maximum adjusted p-value to report. Default `1.0`.

- background_tcrs_df:

  Optional data.frame for background generation. If NULL, uses
  `query_tcrs_df` (which must then include `ja`, `jb`, and nucseq
  columns).

- num_random_samples:

  Integer. Number of random background samples. Default `50000L`.

- fixup_alleles:

  Logical. If TRUE, optimize allele assignments in background TCRs.
  Default `TRUE`.

## Value

A data.frame with columns `tcrdist`, `pvalue_adj`, `fdr_value`,
`query_index` (0-based), `db_index` (0-based), plus query and db TCR
information. Sorted by `pvalue_adj`.

## See also

[`match_tcrs_to_db`](https://shihanli92.github.io/tcrdistR/reference/match_tcrs_to_db.md),
[`strict_single_chain_match_tcrs_to_db`](https://shihanli92.github.io/tcrdistR/reference/strict_single_chain_match_tcrs_to_db.md)

## Examples

``` r
if (FALSE) { # \dontrun{
matches <- find_significant_tcrdist_matches(
    query_tcrs_df = query_df,
    db_tcrs_df = db_df,
    organism = "human",
    adjusted_pvalue_threshold = 0.05
)
} # }
```
