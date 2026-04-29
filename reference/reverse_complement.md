# Reverse complement of a nucleotide sequence

Reverses a nucleotide sequence and replaces each base with its
complement using `BASE_PARTNER`. Supports standard bases (a/c/g/t,
A/C/G/T), IUPAC ambiguity codes, and the gap character (".").

## Usage

``` r
reverse_complement(seq)
```

## Arguments

- seq:

  Character string. The nucleotide sequence to reverse-complement.

## Value

A character string of the same length as `seq`, reversed and
complemented.

## Examples

``` r
reverse_complement("ACGT")   # "ACGT"
#> [1] "ACGT"
reverse_complement("AACGT")  # "ACGTT"
#> [1] "ACGTT"
```
