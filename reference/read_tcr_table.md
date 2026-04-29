# Read a delimited file containing TCR repertoire data

Reads a CSV or TSV file of TCR clonotypes and returns a data.frame with
tcrdistR canonical column names (`va`, `cdr3a`, `vb`, `cdr3b`, and
optionally `ja`, `jb`). Common naming conventions (e.g. tcrdist3/DASH
format) are auto-detected when `col_map` is not provided.

## Usage

``` r
read_tcr_table(file, col_map = NULL, sep = "\t", normalize_genes = TRUE, ...)
```

## Arguments

- file:

  Character string. Path to a delimited text file.

- col_map:

  A named list mapping tcrdistR canonical column names to the column
  names used in `file`. For example,
  `list(va = "v_a_gene", cdr3a = "cdr3_a_aa")`. If `NULL` (the default),
  the function attempts to auto-detect the naming pattern.

- sep:

  Character. Field separator passed to
  [`read.delim`](https://rdrr.io/r/utils/read.table.html). Defaults to
  `"\t"` (tab-separated).

- normalize_genes:

  Logical. If `TRUE` (default), V and J gene names that lack an allele
  suffix (e.g. `TRAV1-1`) are normalised by appending `*01`.

- ...:

  Additional arguments passed to
  [`read.delim`](https://rdrr.io/r/utils/read.table.html).

## Value

A `data.frame` with tcrdistR canonical column names. All gene and CDR3
columns are character vectors (never factors).

## See also

[`read_airr`](https://shihanli92.github.io/tcrdistR/reference/read_airr.md),
[`read_adaptive`](https://shihanli92.github.io/tcrdistR/reference/read_adaptive.md),
[`read_10x`](https://shihanli92.github.io/tcrdistR/reference/read_10x.md),
[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md)

## Examples

``` r
# \donttest{
# Read a tcrdist3-style TSV file
# df <- read_tcr_table("dash_human.tsv")

# Read a CSV with custom column mapping
# df <- read_tcr_table("my_data.csv", sep = ",",
#     col_map = list(va = "V_alpha", cdr3a = "CDR3_alpha",
#                    vb = "V_beta",  cdr3b = "CDR3_beta"))
# }
```
