# Load and cache CD8 logistic regression model parameters

Reads `inst/extdata/cd8_logreg_params_A.txt` and
`inst/extdata/cd8_logreg_params_B.txt`. Each file has 7 header lines
with integer parameters (key-value pairs), followed by tag-weight pairs
for the logistic regression model. The last tag is `BIAS`.

## Usage

``` r
.load_cd8_logreg_models()
```

## Value

A named list with elements `A` and `B`, each containing weights, gene
indexers, and model parameters.
