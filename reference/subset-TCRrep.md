# Subset a TCRrep object by row indices

Subsets a
[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md) by
clonotype indices. Filters `clone_df` and all associated distance
matrices simultaneously. KNN matrices and meta-clonotypes are cleared
since they reference the original indexing.

## Usage

``` r
# S4 method for class 'TCRrep,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]
```

## Arguments

- x:

  A `TCRrep` object.

- i:

  Row indices: logical, integer, or character (rownames).

- j:

  Ignored.

- ...:

  Ignored.

- drop:

  Ignored.

## Value

A new `TCRrep` object with the selected clonotypes.

## Examples

``` r
# \donttest{
data(dash)
rep <- TCRrep(dash, organism = "mouse", compute_distances = TRUE)
#> deduplicate: 1924 -> 1888 clones
pa <- rep[rep@clone_df$epitope == "PA", ]
pa
#> TCRrep object: 321 clonotypes
#>   organism: mouse
#>   chains: AB
#>   metric: tcrdist
#>   distances: paired(dense)
# }
```
