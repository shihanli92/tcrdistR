# Resample shuffled TCR chains

R wrapper for `rcpp_resample_shuffled_tcr_chains`. Generates random TCR
chains by recombining junction segments from existing TCRs.

## Usage

``` r
.resample_shuffled_tcr_chains(
  organism,
  num_samples,
  chain,
  junctions_df,
  preserve_vj_pairings = FALSE,
  max_attempts = 100L * num_samples
)
```

## Arguments

- organism:

  Character string. Organism identifier.

- num_samples:

  Integer. Number of resampled chains to generate.

- chain:

  Character string. Either `"A"` or `"B"`.

- junctions_df:

  A data.frame as produced by `.parse_tcr_junctions`.

- preserve_vj_pairings:

  Logical. If `TRUE`, sample chimeras preserving V-gene and J-gene
  pairings.

## Value

A data.frame with columns `v_gene, j_gene, cdr3, cdr3_nucseq`.
