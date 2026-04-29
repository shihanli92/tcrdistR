# Convert Adaptive Biotechnologies gene names to IMGT format

Transforms gene names from Adaptive's naming convention (e.g.
`TCRBV05-01*01`) to standard IMGT format (e.g. `TRBV5-1*01`). Handles
the prefix change (`TCR` to `TR`), and leading-zero stripping from both
family and sub-family numbers.

## Usage

``` r
.convert_adaptive_gene(gene_names)
```

## Arguments

- gene_names:

  Character vector of Adaptive-format gene names.

## Value

Character vector of the same length with IMGT-formatted names. NAs and
unresolved values are returned as `NA_character_`.
