# List available organisms in the gene database

Returns the names of all organisms available in the bundled gene
database. These are the valid values for the `organism` parameter in
functions like
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)
and
[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md).

## Usage

``` r
list_organisms()
```

## Value

A character vector of organism names, sorted alphabetically.

## See also

[`load_gene_database`](https://shihanli92.github.io/tcrdistR/reference/load_gene_database.md)

## Examples

``` r
list_organisms()
#> [1] "human"     "human_gd"  "human_ig"  "mouse"     "mouse_gd"  "mouse_ig" 
#> [7] "rhesus"    "rhesus_gd"
```
