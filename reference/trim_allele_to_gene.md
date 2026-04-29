# Trim allele designation from a gene ID

Removes the allele suffix (everything from `*` onward) from a TCR or BCR
gene identifier, returning only the gene-level name.

## Usage

``` r
trim_allele_to_gene(gene_id)
```

## Arguments

- gene_id:

  Character string. A gene identifier containing an allele designation,
  e.g. `"TRAV1*01"` or `"IGHV1-18*03"`.

## Value

A character string with the allele suffix stripped, e.g. `"TRAV1"` or
`"IGHV1-18"`.

## See also

[`load_gene_database`](https://shihanli92.github.io/tcrdistR/reference/load_gene_database.md)

## Examples

``` r
trim_allele_to_gene("TRAV1*01")    # "TRAV1"
#> [1] "TRAV1"
trim_allele_to_gene("IGHV1-18*03") # "IGHV1-18"
#> [1] "IGHV1-18"
```
