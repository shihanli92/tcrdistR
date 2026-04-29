# Hamming distance between two CDR3 sequences

Counts the number of positional mismatches between two equal-length
amino acid sequences. Returns `-1L` if the sequences have different
lengths.

## Usage

``` r
hamming_distance(a, b)
```

## Arguments

- a:

  Character string. First CDR3 amino acid sequence.

- b:

  Character string. Second CDR3 amino acid sequence.

## Value

An integer: the number of mismatches, or `-1L` if lengths differ.

## See also

[`hamming_matrix`](https://shihanli92.github.io/tcrdistR/reference/hamming_matrix.md),
[`weighted_cdr3_distance`](https://shihanli92.github.io/tcrdistR/reference/weighted_cdr3_distance.md)

## Examples

``` r
hamming_distance("CASSI", "CASSK")  # 1
#> [1] 1
hamming_distance("CASSI", "CASSI")  # 0
#> [1] 0
hamming_distance("CASSI", "CASSILY")  # -1 (different lengths)
#> [1] -1
```
