# Pairwise Hamming distance matrix for CDR3 sequences

Computes an N x N integer matrix of Hamming distances between CDR3
sequences. For pairs of equal length, the distance is the number of
mismatches. For pairs of unequal length, the distance is set to the
length of the longer sequence (maximum penalty).

## Usage

``` r
hamming_matrix(cdr3_seqs)
```

## Arguments

- cdr3_seqs:

  Character vector. CDR3 amino acid sequences.

## Value

An integer matrix of dimensions N x N.

## See also

[`hamming_distance`](https://shihanli92.github.io/tcrdistR/reference/hamming_distance.md),
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)

## Examples

``` r
seqs <- c("CASSI", "CASSK", "CASRL")
hamming_matrix(seqs)
#>      [,1] [,2] [,3]
#> [1,]    0    1    2
#> [2,]    1    0    2
#> [3,]    2    2    0
```
