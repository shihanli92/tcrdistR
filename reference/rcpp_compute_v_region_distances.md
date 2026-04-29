# Compute pairwise V-region distances (C++ implementation)

Computes all pairwise `blosum_sequence_distance` values for a set of
V-gene CDR loop sequences, weighted by `weight_v_region`. Matches the
inner loop of `compute_all_v_region_distances()` in
`rconga/R/tcr_distances.R`.

## Usage

``` r
rcpp_compute_v_region_distances(
  gene_ids,
  loop_seqs,
  weight_v_region = 1,
  gap_penalty_v_region = 4
)
```

## Arguments

- gene_ids:

  Character vector of length N. V-gene allele identifiers (e.g.
  `"TRAV1-1*01"`). Used as row/column names of the result.

- loop_seqs:

  Character vector of length N. Pre-aligned merged CDR loop sequences
  corresponding to each gene in `gene_ids`. All sequences must have the
  same character length.

- weight_v_region:

  Numeric. Multiplied by each raw `blosum_sequence_distance` value.
  Default 1.0 matches `WEIGHT_V_REGION`.

- gap_penalty_v_region:

  Numeric. Forwarded to `rcpp_blosum_sequence_distance` as
  `gap_penalty`. Default 4.0 matches `GAP_PENALTY_V_REGION`.

## Value

A named symmetric N x N `NumericMatrix`.
`result[i, j] = weight_v_region * blosum_sequence_distance(loop_seqs[i], loop_seqs[j])`.

## Details

The sequences in `loop_seqs` are pre-aligned merged CDR loop strings
(CDRs joined by `" "`), as produced by the R gene-database loader.
Spaces mark position padding and are skipped (both positions must be
space).

## Examples

``` r
if (FALSE) { # \dontrun{
  ids  <- c("TRAV1-1*01", "TRAV1-2*01")
  seqs <- c("ACDE FGHI", "ACDE FGHK")
  rcpp_compute_v_region_distances(ids, seqs)
} # }
```
