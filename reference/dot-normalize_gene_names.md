# Normalize V/J gene names by appending default allele suffix

If a gene name lacks a `*` allele suffix, appends `*01`. NAs and empty
strings are passed through unchanged.

## Usage

``` r
.normalize_gene_names(gene_names)
```

## Arguments

- gene_names:

  Character vector of gene names.

## Value

Character vector of the same length with allele suffixes ensured.
