# TCR repertoire clonality

Computes clonality as 1 minus the normalized Shannon entropy:
`1 - H / log(S)` where `H = -sum(p * log(p))` and `S` is the number of
unique clonotypes. Values near 0 indicate a uniform (diverse)
repertoire; values near 1 indicate a dominated repertoire.

## Usage

``` r
tcr_clonality(counts)
```

## Arguments

- counts:

  Integer vector. Clonotype counts (positive integers).

## Value

A numeric scalar between 0 and 1.

## See also

[`tcr_diversity`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md),
[`tcr_richness`](https://shihanli92.github.io/tcrdistR/reference/tcr_richness.md),
[`tcr_shannon_entropy`](https://shihanli92.github.io/tcrdistR/reference/tcr_shannon_entropy.md)

## Examples

``` r
# Uniform: clonality near 0
tcr_clonality(rep(10, 5))
#> [1] -2.220446e-16

# Dominated: clonality near 1
tcr_clonality(c(1000, 1, 1))
#> [1] 0.985631
```
