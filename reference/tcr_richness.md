# TCR repertoire richness

Returns the number of unique clonotypes (species richness).

## Usage

``` r
tcr_richness(counts)
```

## Arguments

- counts:

  Integer vector. Clonotype counts (positive integers).

## Value

An integer: the number of unique clonotypes.

## See also

[`tcr_diversity`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md),
[`tcr_clonality`](https://shihanli92.github.io/tcrdistR/reference/tcr_clonality.md)

## Examples

``` r
tcr_richness(c(10, 5, 3, 1, 1))  # 5
#> [1] 5
```
