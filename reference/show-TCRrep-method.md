# Show method for TCRrep objects

Prints a compact summary of a
[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md)
object, including the number of clonotypes, organism, chain combination,
distance metric, distance computation status, and k-nearest-neighbour
status when available.

## Usage

``` r
# S4 method for class 'TCRrep'
show(object)
```

## Arguments

- object:

  A `TCRrep` object.

## Value

Returns `invisible(object)`.

## Examples

``` r
# \donttest{
tcrs <- data.frame(
    va    = c("TRAV1-1*01", "TRAV1-1*01"),
    cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
    vb    = c("TRBV19*01", "TRBV19*01"),
    cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
    stringsAsFactors = FALSE
)
obj <- TCRrep(tcrs, organism = "human")
show(obj)
#> TCRrep object: 2 clonotypes
#>   organism: human
#>   chains: AB
#>   metric: tcrdist
#>   distances: not computed
# }
```
