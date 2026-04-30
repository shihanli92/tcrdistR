# Identify quasi-public meta-clonotypes

Finds TCR sequences that are shared across multiple subjects
(individuals) based on TCRdist neighborhoods. A meta-clonotype is
defined by a center TCR and a radius: all TCRs within the radius are
considered part of the meta-clonotype.

## Usage

``` r
find_meta_clonotypes(
  tcr_df,
  organism,
  radius = NULL,
  background_df = NULL,
  ctrl_bkgd = 1e-05,
  max_radius = 50L,
  min_nsubject = 2L,
  subject_col = "subject",
  dist_matrix = NULL
)
```

## Arguments

- tcr_df:

  Data.frame with TCR columns (`va`, `vb`, `cdr3a`, `cdr3b`) and a
  subject column.

- organism:

  Character string (`"human"` or `"mouse"`).

- radius:

  Numeric or integer vector. TCRdist radius for neighborhood definition.
  If a single value, used for all TCRs. If `NULL`, per-TCR radii are
  computed from background ECDF.

- background_df:

  Data.frame. Background TCRs for radius computation. Required if
  `radius` is `NULL`.

- ctrl_bkgd:

  Numeric. Background proportion threshold for automatic radius. Default
  `1e-5`.

- max_radius:

  Numeric. Maximum radius to consider. Default `50L`.

- min_nsubject:

  Integer. Minimum number of distinct subjects in a neighborhood.
  Default `2L`.

- subject_col:

  Character string. Column name for subject IDs. Default `"subject"`.

- dist_matrix:

  Optional precomputed distance matrix. If provided, `tcr_df` and
  `organism` are not used for distance computation.

## Value

A data.frame with one row per meta-clonotype center:

- `center_index`:

  Row index in `tcr_df`.

- `va, cdr3a, vb, cdr3b`:

  Center TCR sequence.

- `radius`:

  TCRdist radius used.

- `K_neighbors`:

  Total neighbors within radius.

- `nsubject`:

  Number of distinct subjects.

- `neighbor_indices`:

  Comma-separated neighbor row indices.

## See also

[`summarize_meta_clonotype`](https://shihanli92.github.io/tcrdistR/reference/summarize_meta_clonotype.md),
[`find_clumping`](https://shihanli92.github.io/tcrdistR/reference/find_clumping.md)

## Examples

``` r
if (FALSE) { # \dontrun{
meta <- find_meta_clonotypes(tcr_df, "human", radius = 30,
                              subject_col = "subject")
} # }
```
