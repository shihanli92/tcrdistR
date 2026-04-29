# Read 10X Genomics VDJ contig annotations

Reads a 10X Genomics cellranger VDJ output file (typically
`filtered_contig_annotations.csv`) and returns a paired alpha-beta TCR
data frame suitable for
[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md).
Each row in the 10X format represents a single contig; this function
filters, deduplicates, and pairs alpha and beta chains by cell barcode.

## Usage

``` r
read_10x(
  file,
  pair_by = "barcode",
  productive_only = TRUE,
  normalize_genes = TRUE
)
```

## Arguments

- file:

  Character string. Path to a 10X VDJ CSV file (e.g.
  `filtered_contig_annotations.csv`).

- pair_by:

  Character string. Column name used to pair TRA and TRB contigs from
  the same cell. Defaults to `"barcode"`.

- productive_only:

  Logical. If `TRUE` (default), only productive contigs are retained
  (`productive == "True"`).

- normalize_genes:

  Logical. If `TRUE` (default), gene names lacking an allele suffix
  (e.g. `"TRAV1-1"`) are appended with `*01`.

## Value

A `data.frame` with columns `barcode`, `va`, `ja`, `cdr3a`, `vb`, `jb`,
`cdr3b`. Only cells with both a TRA and TRB contig are included (inner
join).

## Details

When a cell has multiple contigs for the same chain, the contig with the
highest UMI count is retained (if the `umis` column is present);
otherwise the first contig is kept.

## See also

[`read_tcr_table`](https://shihanli92.github.io/tcrdistR/reference/read_tcr_table.md),
[`read_airr`](https://shihanli92.github.io/tcrdistR/reference/read_airr.md),
[`read_adaptive`](https://shihanli92.github.io/tcrdistR/reference/read_adaptive.md)

## Examples

``` r
# \donttest{
# df <- read_10x("filtered_contig_annotations.csv")
# obj <- TCRrep(df, organism = "human")
# }
```
