# Generalized Simpson's diversity (Hill numbers framework)

Computes the order-`r` generalized Simpson's entropy for a vector of
clonotype counts. For `order=2`, this is the classical Simpson's
diversity index. Higher orders give more weight to dominant clonotypes.

## Usage

``` r
tcr_diversity(counts, order = 2L, ci = TRUE, alpha = 0.05)
```

## Arguments

- counts:

  Integer vector. Clonotype counts (positive integers).

- order:

  Integer. Order of the diversity index. Default `2L` (standard
  Simpson's).

- ci:

  Logical. If `TRUE` (default), compute confidence interval.

- alpha:

  Numeric. Significance level for CI. Default `0.05`.

## Value

A named list:

- `entropy`:

  Numeric. The diversity index Z_r, between 0 and 1.

- `effective_number`:

  Numeric. Hill number (effective species).

- `ci_lower`:

  Numeric. Lower CI bound (if `ci=TRUE`).

- `ci_upper`:

  Numeric. Upper CI bound (if `ci=TRUE`).

- `order`:

  Integer. The order used.

## Details

The diversity index \\Z_r\\ is the probability that \\r\\ randomly drawn
individuals (without replacement) all belong to the same species:

\$\$Z_r = \sum_i \frac{c_i}{n} \prod\_{k=1}^{r-1} \frac{c_i - k}{n -
k}\$\$

where \\c_i\\ is the count of clonotype \\i\\ and \\n = \sum c_i\\. The
returned entropy is \\1 - Z_r\\ (probability of drawing \\r\\ different
species). For \\r=2\\ this is the classical Simpson's diversity.

The effective number of species (Hill number) is: \\D = 1 /
Z_r^{1/(r-1)}\\.

Confidence intervals use the delta method with the multinomial variance
of \\Z_r\\. For \\r=2\\, the exact gradient \\\partial Z / \partial p_i
= 2 p_i\\ is used; for higher orders, an approximation \\r \cdot
p_i^{r-1}\\ is applied.

## References

Hill, M. O. (1973). Diversity and evenness: a unifying notation and its
consequences. *Ecology*, 54(2), 427–432.

Jost, L. (2006). Entropy and diversity. *Oikos*, 113(2), 363–375.

## See also

[`tcr_fuzzy_diversity`](https://shihanli92.github.io/tcrdistR/reference/tcr_fuzzy_diversity.md),
[`tcr_richness`](https://shihanli92.github.io/tcrdistR/reference/tcr_richness.md),
[`tcr_clonality`](https://shihanli92.github.io/tcrdistR/reference/tcr_clonality.md)

## Examples

``` r
# Uniform distribution: maximum diversity
tcr_diversity(rep(10, 5))
#> $entropy
#> [1] 0.8163265
#> 
#> $effective_number
#> [1] 5.444444
#> 
#> $order
#> [1] 2
#> 
#> $ci_lower
#> [1] 0.8163265
#> 
#> $ci_upper
#> [1] 0.8163265
#> 

# Single dominant clonotype: low diversity
tcr_diversity(c(100, 1, 1, 1))
#> $entropy
#> [1] 0.05768132
#> 
#> $effective_number
#> [1] 1.061212
#> 
#> $order
#> [1] 2
#> 
#> $ci_lower
#> [1] 0
#> 
#> $ci_upper
#> [1] 0.1201095
#> 
```
