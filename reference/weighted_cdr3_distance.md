# Weighted CDR3 distance between two CDR3 sequences

Computes the TCRdist CDR3 component distance between two amino acid CDR3
sequences. Accounts for length differences via gap penalties and uses
BSD4 substitution scores on the aligned positions. Dispatches to a C++
implementation for performance.

## Usage

``` r
weighted_cdr3_distance(
  seq1,
  seq2,
  weight = WEIGHT_CDR3_REGION,
  gap_penalty = GAP_PENALTY_CDR3_REGION
)
```

## Arguments

- seq1:

  Character string of length 1. First CDR3 amino acid sequence. Must be
  non-empty.

- seq2:

  Character string of length 1. Second CDR3 amino acid sequence. Must be
  non-empty.

- weight:

  Integer. Weight applied to the alignment score component. Defaults to
  `WEIGHT_CDR3_REGION` (3L).

- gap_penalty:

  Integer. Penalty per gap character in length differences. Defaults to
  `GAP_PENALTY_CDR3_REGION` (12L).

## Value

A numeric scalar: the weighted CDR3 distance.

## See also

[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md),
[`bsd4_matrix`](https://shihanli92.github.io/tcrdistR/reference/bsd4_matrix.md)

## Examples

``` r
# \donttest{
weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSYEQYF")
#> [1] 12
weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSSYEQYF")  # 0
#> [1] 0
# }
```
