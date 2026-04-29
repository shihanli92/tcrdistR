# Retrieve the BSD4 substitution matrix

Returns the BSD4 (BLOSUM-derived substitution distance 4) matrix used
internally for amino acid distance scoring in TCRdist calculations. The
matrix is built in C++ and returned as a named numeric matrix in R.

## Usage

``` r
bsd4_matrix()
```

## Value

A named numeric matrix of amino acid substitution distances. Row and
column names are the 20 standard single-letter amino acid codes in
`AMINO_ACIDS` order.

## See also

[`weighted_cdr3_distance`](https://shihanli92.github.io/tcrdistR/reference/weighted_cdr3_distance.md),
[`AMINO_ACIDS`](https://shihanli92.github.io/tcrdistR/reference/AMINO_ACIDS.md)

## Examples

``` r
# \donttest{
bsd4 <- bsd4_matrix()
bsd4["A", "A"]  # 0
#> [1] 0
bsd4["A", "G"]  # small positive value
#> [1] 4
# }
```
