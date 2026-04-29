# Find alternate alleles for a batch of TCRs

R wrapper for `rcpp_find_alternate_alleles_batch`. For each TCR, tries
alternate alleles. If an alternate allele is consistently better (count
\>= `min_better_count` AND ratio \>= `min_better_ratio`:1), swaps it
into the data.

## Usage

``` r
.find_alternate_alleles_for_tcrs(
  organism,
  tcr_df,
  min_better_ratio = 10,
  min_better_count = 5L,
  min_improvement = 2L,
  verbose = FALSE
)
```

## Arguments

- organism:

  Character string. Organism identifier.

- tcr_df:

  A data.frame with columns
  `va, ja, cdr3a_nucseq, vb, jb, cdr3b_nucseq`.

- min_better_ratio:

  Numeric. Minimum ratio for the alternate allele to be allowed.

- min_better_count:

  Integer. Minimum count for the alternate allele to be allowed.

- min_improvement:

  Integer. Minimum improvement in match count.

- verbose:

  Logical. If `TRUE`, print progress information.

## Value

A data.frame of the same structure as `tcr_df` with potentially updated
gene names.
