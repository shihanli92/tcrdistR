# Gini index of a clonotype distribution

Computes the Gini coefficient measuring inequality in clonotype
abundances. Values near 0 indicate a uniform (equal) distribution;
values near 1 indicate extreme inequality (one dominant clone).

## Usage

``` r
tcr_gini(counts)
```

## Arguments

- counts:

  Integer vector. Clonotype counts (positive integers).

## Value

A numeric scalar between 0 and 1.

## See also

[`tcr_shannon_entropy`](https://shihanli92.github.io/tcrdistR/reference/tcr_shannon_entropy.md),
[`tcr_clonality`](https://shihanli92.github.io/tcrdistR/reference/tcr_clonality.md),
[`tcr_diversity`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md)

## Examples

``` r
# Uniform: Gini near 0
tcr_gini(rep(10, 5))
#> [1] 0

# Dominated: Gini near 1
tcr_gini(c(10000, 1, 1))
#> [1] 0.6664667
```
