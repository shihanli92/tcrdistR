# Parse junction regions for a batch of TCRs

Batch wrapper around `.analyze_junction`. Processes each TCR and returns
a data.frame with junction analysis results.

## Usage

``` r
.parse_tcr_junctions(organism, tcr_df)
```

## Arguments

- organism:

  Character string. Organism identifier.

- tcr_df:

  A data.frame with columns
  `va, ja, cdr3a, cdr3a_nucseq, vb, jb, cdr3b, cdr3b_nucseq`.

## Value

A data.frame with junction info columns needed for resampling.
