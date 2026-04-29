# Read AIRR Rearrangement TSV format

Reads a file in the AIRR Community Standard Rearrangement TSV format and
returns a paired alpha-beta TCR data frame suitable for
[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md).
Each row in the AIRR format represents a single rearrangement; this
function pairs alpha and beta chains by a shared identifier (typically
`cell_id`).

## Usage

``` r
read_airr(
  file,
  pair_by = "cell_id",
  productive_only = TRUE,
  normalize_genes = TRUE
)
```

## Arguments

- file:

  Character string. Path to an AIRR Rearrangement TSV file.

- pair_by:

  Character string. Column name used to pair TRA and TRB rearrangements
  from the same cell. Defaults to `"cell_id"`.

- productive_only:

  Logical. If `TRUE` (default), only productive rearrangements are
  retained. The `productive` column may contain logical `TRUE` or
  character `"T"`/`"TRUE"`.

- normalize_genes:

  Logical. If `TRUE` (default), gene names lacking an allele suffix
  (e.g. `"TRAV1-1"`) are appended with `*01`.

## Value

A `data.frame` with columns `cell_id`, `va`, `ja`, `cdr3a`, `vb`, `jb`,
`cdr3b`. Only cells with both a TRA and TRB rearrangement are included
(inner join).

## Details

Multi-value gene calls (comma-separated, e.g. `"TRAV1-1*01,TRAV1-2*01"`)
are resolved by taking the first listed allele.

## See also

[`read_tcr_table`](https://shihanli92.github.io/tcrdistR/reference/read_tcr_table.md),
[`read_adaptive`](https://shihanli92.github.io/tcrdistR/reference/read_adaptive.md),
[`read_10x`](https://shihanli92.github.io/tcrdistR/reference/read_10x.md)

## Examples

``` r
# \donttest{
# df <- read_airr("rearrangements.tsv")
# obj <- TCRrep(df, organism = "human")
# }
```
