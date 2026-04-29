# Count nucleotide prefix matches between two sequences (C++ implementation)

Linear scan prefix matching: +1 per match, mismatch_score per mismatch.
Returns the length of the best-scoring prefix (the position at which the
cumulative score was highest).

## Usage

``` r
rcpp_count_nuc_matches(a, b, mismatch_score = -4L)
```

## Arguments

- a:

  Character string. First nucleotide sequence.

- b:

  Character string. Second nucleotide sequence.

- mismatch_score:

  Integer. Score added per mismatch position. Default -4.

## Value

Integer: length of best-scoring prefix (0 if either string is empty).

## Examples

``` r
if (FALSE) { # \dontrun{
  rcpp_count_nuc_matches("ATCGATCG", "ATCGTTCG")
} # }
```
