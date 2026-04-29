# Extract edge list from a sparse distance matrix

Converts upper triangle of a symmetric sparse matrix to an edge list
with distance values.

## Usage

``` r
.sparse_to_edgelist(sp)
```

## Arguments

- sp:

  A `dgCMatrix` (symmetric sparse distance matrix).

## Value

A `data.frame` with columns `from`, `to`, `distance` (all entries from
upper triangle with dist \> 0).
