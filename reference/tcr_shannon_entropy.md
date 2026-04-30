# Shannon entropy of a clonotype distribution

Computes Shannon entropy `H = -sum(p * log(p, base))` for a vector of
clonotype counts. Natural log (nats) by default; set `base = 2` for
bits.

## Usage

``` r
tcr_shannon_entropy(counts, base = exp(1))
```

## Arguments

- counts:

  Integer vector. Clonotype counts (positive integers).

- base:

  Numeric. Logarithm base. Default `exp(1)` (natural log). Use `2` for
  bits.

## Value

A numeric scalar \>= 0.

## See also

[`tcr_clonality`](https://shihanli92.github.io/tcrdistR/reference/tcr_clonality.md),
[`tcr_diversity`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md),
[`tcr_gini`](https://shihanli92.github.io/tcrdistR/reference/tcr_gini.md)

## Examples

``` r
# Uniform: H = log(S)
tcr_shannon_entropy(rep(10, 5))
#> [1] 1.609438

# Single species: H = 0
tcr_shannon_entropy(c(100))
#> [1] 0
```
