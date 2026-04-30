# Repertoire overlap metrics

Computes overlap metrics between two named count vectors representing
clonotype abundances from two samples. Available metrics:

## Usage

``` r
tcr_repertoire_overlap(
  counts_a,
  counts_b,
  metrics = c("jaccard", "morisita_horn", "overlap_coef")
)
```

## Arguments

- counts_a, counts_b:

  Named numeric vectors of clonotype counts. Names are clonotype
  identifiers.

- metrics:

  Character vector. Which metrics to compute. Default computes all
  three.

## Value

A named list with the requested metrics, each a scalar in `[0, 1]`.

## Details

- Jaccard index:

  `|A intersect B| / |A union B|` — proportion of shared clonotype
  species (presence/absence).

- Morisita-Horn:

  `2 * sum(p_a * p_b) / (sum(p_a^2) + sum(p_b^2))` — abundance-weighted
  overlap.

- Overlap coefficient:

  `|A intersect B| / min(|A|, |B|)` — overlap normalized by the smaller
  set.

## See also

[`tcr_diversity`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md),
[`tcr_clonality`](https://shihanli92.github.io/tcrdistR/reference/tcr_clonality.md)

## Examples

``` r
a <- c(clone1 = 10, clone2 = 5, clone3 = 1)
b <- c(clone1 = 8, clone3 = 3, clone4 = 2)
tcr_repertoire_overlap(a, b)
#> $jaccard
#> [1] 0.5
#> 
#> $morisita_horn
#> [1] 0.8420231
#> 
#> $overlap_coef
#> [1] 0.6666667
#> 
```
