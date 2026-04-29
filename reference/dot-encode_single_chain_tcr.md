# Encode a single TCR chain as a feature vector for CD8 logistic regression

Creates a one-hot/k-hot feature vector encoding V-gene, J-gene, CDR3
length, and CDR3 amino acid composition at N-terminal, C-terminal, and
middle positions.

## Usage

``` r
.encode_single_chain_tcr(vgene, jgene, cdr3, model_params)
```

## Arguments

- vgene:

  Character. V-gene name with allele.

- jgene:

  Character. J-gene name with allele.

- cdr3:

  Character. CDR3 amino acid sequence.

- model_params:

  List. Model parameters from `.load_cd8_logreg_models`.

## Value

Numeric vector of length NTOT + 1 (features + bias).
