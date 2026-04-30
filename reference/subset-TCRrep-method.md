# Subset a TCRrep using non-standard evaluation

Filter clonotypes using column expressions evaluated against `clone_df`.
This is the recommended way to filter a TCRrep:

## Usage

``` r
# S4 method for class 'TCRrep'
subset(x, subset, ...)
```

## Arguments

- x:

  A `TCRrep` object.

- subset:

  A logical expression evaluated in the context of `x@clone_df`. Column
  names can be used directly.

- ...:

  Ignored.

## Value

A new `TCRrep` object with only matching clonotypes.

## Details

    rep |> subset(epitope == "PA")

## Examples

``` r
# \donttest{
data(dash)
rep <- TCRrep(dash, organism = "mouse", compute_distances = TRUE)
#> deduplicate: 1924 -> 1888 clones
pa <- subset(rep, epitope == "PA")
pa
#> TCRrep object: 321 clonotypes
#>   organism: mouse
#>   chains: AB
#>   metric: tcrdist
#>   distances: paired(dense)
# }
```
