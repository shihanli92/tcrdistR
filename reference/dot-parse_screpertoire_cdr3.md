# Parse scRepertoire CTaa/CTnt column into alpha/beta CDR3 columns

Parses scRepertoire's concatenated CDR3 format: `cdr3a_cdr3b`
(underscore-separated). Multi-chain cells (semicolon-separated) use the
first chain per locus.

## Usage

``` r
.parse_screpertoire_cdr3(cdr3_col)
```

## Arguments

- cdr3_col:

  Character vector. The `CTaa` or `CTnt` column.

## Value

A `data.frame` with columns `alpha` and `beta`.
