# Test association between a variable and TCR neighborhoods

For each TCR, tests whether a categorical variable is non-randomly
distributed among its TCRdist neighbors compared to the full repertoire.
Supports Fisher's exact test (binary variables) and chi-squared test
(multi-category variables).

## Usage

``` r
neighborhood_test(
  tcr_df,
  organism,
  variable,
  radius = 50,
  test = c("fisher", "chisq"),
  p_adjust_method = "BH"
)
```

## Arguments

- tcr_df:

  Data.frame with TCR columns.

- organism:

  Character string (`"human"` or `"mouse"`).

- variable:

  Character or factor vector of length `nrow(tcr_df)`. The categorical
  variable to test.

- radius:

  Numeric. Maximum TCRdist for neighborhood membership. Default `50`.

- test:

  Character string. `"fisher"` (default, for binary) or `"chisq"` (for
  multi-category).

- p_adjust_method:

  Character string. Method for
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"BH"` (Benjamini-Hochberg).

## Value

A data.frame with one row per TCR and columns:

- `index`:

  Row index in `tcr_df`.

- `n_neighbors`:

  Number of neighbors within radius.

- `p_value`:

  Raw test p-value.

- `p_adjusted`:

  Adjusted p-value.

- `odds_ratio`:

  Odds ratio (Fisher only, NA for chi-sq).

## See also

[`tcrdist_hclust`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_hclust.md),
[`cluster_tcrs`](https://shihanli92.github.io/tcrdistR/reference/cluster_tcrs.md),
[`find_clumping`](https://shihanli92.github.io/tcrdistR/reference/find_clumping.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- neighborhood_test(tcr_df, "human",
                             variable = tcr_df$epitope, radius = 50)
significant <- result[result$p_adjusted < 0.05, ]
} # }
```
