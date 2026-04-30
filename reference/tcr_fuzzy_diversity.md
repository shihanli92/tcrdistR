# TCR-aware fuzzy diversity

Computes diversity accounting for sequence similarity: two clonotypes
are considered "the same" if their TCRdist is within `threshold`. This
gives lower diversity for repertoires with many similar sequences.

## Usage

``` r
tcr_fuzzy_diversity(
  tcr_df,
  organism,
  threshold = 50,
  order = 2L,
  counts = NULL
)
```

## Arguments

- tcr_df:

  Data.frame with TCR columns (`va`, `vb`, `cdr3a`, `cdr3b`).

- organism:

  Character string (`"human"` or `"mouse"`).

- threshold:

  Numeric. Distance threshold for considering two TCRs as similar.
  Default `50`.

- order:

  Integer. Diversity order. Default `2L`.

- counts:

  Integer vector. Clonotype counts (one per row of `tcr_df`). If `NULL`
  (default), all counts are 1.

## Value

A named list:

- `fuzzy_diversity`:

  Numeric. Fuzzy diversity, between 0 and 1.

- `standard_diversity`:

  Numeric. Standard Simpson's diversity for comparison.

## Details

For `order=2`, the fuzzy Simpson's index is computed analytically:

\$\$Z\_{fuzzy} = \frac{\sum\_{i,j} c_i \cdot c_j \cdot I(d(i,j) \le
threshold)}{(\sum_i c_i)^2}\$\$

where \\I(\cdot)\\ is the indicator function and \\d(i,j)\\ is the
TCRdist between clonotypes \\i\\ and \\j\\. The fuzzy diversity is \\1 -
Z\_{fuzzy}\\. This is always \\\le\\ the standard Simpson's diversity
because merging similar clonotypes increases the concentration.

For higher orders, a sampling-based approximation is used (10,000
draws).

## See also

[`tcr_diversity`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md),
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tcr_fuzzy_diversity(tcr_df, "human", threshold = 50)
} # }
```
