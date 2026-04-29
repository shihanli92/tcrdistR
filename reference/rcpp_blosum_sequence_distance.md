# Aligned sequence distance using BSD4 substitution costs (C++ implementation)

Computes the total BSD4 distance between two pre-aligned amino acid
strings of equal length. Matches `blosum_sequence_distance()` in
`rconga/R/tcr_distances.R` exactly.

## Usage

``` r
rcpp_blosum_sequence_distance(seq1, seq2, gap_penalty = 4)
```

## Arguments

- seq1:

  Character string. Aligned sequence (may contain `" "`, `"."`, or
  `"*"`).

- seq2:

  Character string. Aligned sequence of the same length as `seq1`.

- gap_penalty:

  Numeric. Applied when exactly one position is a gap or stop-codon
  character. Default 4.0 matches `GAP_PENALTY_V_REGION`.

## Value

Numeric scalar: sum of per-position BSD4 distances.

## Details

Position handling:

- Both characters are `" "` (space): position skipped (padding). The
  function stops with an error if only one character is a space.

- Both characters are `"."` (gap): distance 0.

- Both characters are `"*"` (stop codon): distance 0.

- One character is `"."` or `"*"`: `gap_penalty`.

- Both are valid amino acids: BSD4 lookup from the constexpr table.

## Examples

``` r
if (FALSE) { # \dontrun{
  rcpp_blosum_sequence_distance("ACDE", "ACDF")
} # }
```
