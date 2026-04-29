# Read Adaptive Biotechnologies ImmunoSeq TSV file

Imports a beta-chain TCR repertoire from an Adaptive ImmunoSeq export
file (TSV). Gene names are converted from Adaptive naming to IMGT format
and optionally normalized to include an allele suffix.

## Usage

``` r
read_adaptive(file, normalize_genes = TRUE)
```

## Arguments

- file:

  Character string. Path to an Adaptive ImmunoSeq TSV file.

- normalize_genes:

  Logical. If `TRUE` (the default), gene names that lack an allele
  suffix (`*XX`) are appended with `*01` after Adaptive-to-IMGT
  conversion.

## Value

A `data.frame` with columns `vb`, `jb`, `cdr3b`, and optionally `count`.
Gene names are in IMGT format. Rows with empty or NA CDR3 sequences are
removed.

## Details

Two column naming conventions are supported automatically:

- Convention 1 (newer):

  `vGeneName`, `jGeneName`, `aminoAcid`, `count`

- Convention 2 (older):

  `v_resolved` (or `vFamilyName`), `j_resolved` (or `jFamilyName`),
  `aminoAcid`, `count` / `templates` / `reads`

## See also

[`read_tcr_table`](https://shihanli92.github.io/tcrdistR/reference/read_tcr_table.md),
[`read_airr`](https://shihanli92.github.io/tcrdistR/reference/read_airr.md),
[`read_10x`](https://shihanli92.github.io/tcrdistR/reference/read_10x.md)

## Examples

``` r
# \donttest{
# df <- read_adaptive("sample_immunoseq.tsv")
# head(df)
# }
```
