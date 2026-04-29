# Strict CDR3 sequence matching against a database

Finds database entries where the CDR3alpha or CDR3beta amino acid
sequence exactly matches a query TCR.

## Usage

``` r
strict_single_chain_match_tcrs_to_db(tcr_df, organism, db_tcrs_tsvfile = NULL)
```

## Arguments

- tcr_df:

  A data.frame with columns `cdr3a` and `cdr3b`.

- organism:

  Character string (e.g. `"human"`, `"mouse"`).

- db_tcrs_tsvfile:

  Path to a TSV file. If NULL, uses the built-in organism-specific
  single-chain database.

## Value

A list with elements `alpha_matches` and `beta_matches`, each a
data.frame of matching database rows.

## See also

[`match_tcrs_to_db`](https://shihanli92.github.io/tcrdistR/reference/match_tcrs_to_db.md),
[`find_significant_tcrdist_matches`](https://shihanli92.github.io/tcrdistR/reference/find_significant_tcrdist_matches.md)

## Examples

``` r
if (FALSE) { # \dontrun{
hits <- strict_single_chain_match_tcrs_to_db(tcr_df, "human")
hits$alpha_matches
hits$beta_matches
} # }
```
