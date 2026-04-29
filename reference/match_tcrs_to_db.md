# Match TCRs to a literature database

Wrapper around
[`find_significant_tcrdist_matches`](https://shihanli92.github.io/tcrdistR/reference/find_significant_tcrdist_matches.md)
that loads a default or custom database file.

## Usage

``` r
match_tcrs_to_db(
  tcr_df,
  organism,
  db_tcrs_tsvfile = NULL,
  adjusted_pvalue_threshold = 1,
  background_tcrs_df = NULL,
  num_random_samples = 50000L
)
```

## Arguments

- tcr_df:

  A data.frame with columns
  `va, ja, cdr3a, cdr3a_nucseq, vb, jb, cdr3b, cdr3b_nucseq`.

- organism:

  Character string.

- db_tcrs_tsvfile:

  Path to a TSV file with database TCRs. If NULL, uses the built-in
  `new_paired_tcr_db_for_matching_nr.tsv` (human only).

- adjusted_pvalue_threshold:

  Numeric. Default `1.0`.

- background_tcrs_df:

  Optional background data.frame.

- num_random_samples:

  Integer. Default `50000L`.

## Value

A data.frame of significant matches (see
[`find_significant_tcrdist_matches`](https://shihanli92.github.io/tcrdistR/reference/find_significant_tcrdist_matches.md)).

## See also

[`find_significant_tcrdist_matches`](https://shihanli92.github.io/tcrdistR/reference/find_significant_tcrdist_matches.md),
[`strict_single_chain_match_tcrs_to_db`](https://shihanli92.github.io/tcrdistR/reference/strict_single_chain_match_tcrs_to_db.md)

## Examples

``` r
if (FALSE) { # \dontrun{
results <- match_tcrs_to_db(tcr_df, organism = "human")
results[results$pvalue_adj < 0.05, ]
} # }
```
