# Join two TCR datasets by TCRdist threshold

Performs a fuzzy merge of two TCR data.frames: for each row in
`left_df`, finds rows in `right_df` within `radius` TCRdist units.
Similar to a SQL JOIN but using distance instead of exact matching.

## Usage

``` r
tcrdist_join(
  left_df,
  right_df,
  organism,
  radius,
  max_n = 5L,
  type = c("inner", "left"),
  suffix = c("_x", "_y")
)
```

## Arguments

- left_df:

  Data.frame with TCR columns (`va`, `vb`, `cdr3a`, `cdr3b`) plus any
  extra columns.

- right_df:

  Data.frame with TCR columns (`va`, `vb`, `cdr3a`, `cdr3b`) plus any
  extra columns.

- organism:

  Character string (`"human"` or `"mouse"`).

- radius:

  Numeric. Maximum TCRdist for a match.

- max_n:

  Integer. Maximum number of matches per left row. Default `5L`. Closest
  matches are kept.

- type:

  Character string. Join type: `"inner"` (default) keeps only matched
  pairs, `"left"` keeps all left rows (NAs for unmatched).

- suffix:

  Character vector of length 2. Suffixes for disambiguating column
  names. Default `c("_x", "_y")`.

## Value

A data.frame with columns from both sides (suffixed if overlapping) plus
a `tcrdist` column.

## See also

[`tcrdist_rect`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_rect.md),
[`tcrdist_radius_neighbors`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_radius_neighbors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
matches <- tcrdist_join(query_tcrs, reference_tcrs, "human", radius = 50)
} # }
```
