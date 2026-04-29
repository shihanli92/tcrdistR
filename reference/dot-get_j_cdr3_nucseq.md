# Get J-gene CDR3 region nucleotide sequence

Returns the J-gene CDR3 region nucleotide sequence up to (but not
including) the GXG motif.

## Usage

``` r
.get_j_cdr3_nucseq(organism, j_gene)
```

## Arguments

- organism:

  Character string. Organism identifier.

- j_gene:

  Character string. J-gene identifier including allele.

## Value

Character string. Lowercase nucleotide sequence. Returns `""` if the
gene is not found.
