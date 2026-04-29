# Weighted CDR3 distance (C++ implementation)

Direct C++ port of the CDR3 distance algorithm from
`rconga/src/tcrdist_rcpp.cpp`. Uses hard-coded `ntrim = 3`, `ctrim = 2`
(`TRIM_CDR3S = TRUE`) and the fixed gap position formula
`min(6, 3 + (lenshort - 5) FALSE`).

## Usage

``` r
rcpp_weighted_cdr3_distance(
  seq1,
  seq2,
  weight_cdr3_region = 3L,
  gap_penalty_cdr3_region = 12L
)
```

## Arguments

- seq1:

  Character string. First CDR3 amino acid sequence.

- seq2:

  Character string. Second CDR3 amino acid sequence.

- weight_cdr3_region:

  Integer. Multiplied by the raw alignment distance. Default 3 matches
  `WEIGHT_CDR3_REGION`.

- gap_penalty_cdr3_region:

  Integer. Per-residue length-difference penalty. Default 12 matches
  `GAP_PENALTY_CDR3_REGION`.

## Value

Numeric scalar: the weighted CDR3 distance.

## Details

BSD4 values are compiled in as `constexpr` (no matrix argument needed).
The result is numerically identical to
`rconga::rcpp_weighted_cdr3_distance()` for valid inputs.

## Examples

``` r
if (FALSE) { # \dontrun{
  rcpp_weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSYEQYF")
} # }
```
