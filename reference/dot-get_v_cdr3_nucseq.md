# Get V-gene CDR3 region nucleotide sequence

Returns the V-gene CDR3 region nucleotide sequence (lowercase), starting
from the conserved C codon.

## Usage

``` r
.get_v_cdr3_nucseq(organism, v_gene)
```

## Arguments

- organism:

  Character string. Organism identifier.

- v_gene:

  Character string. V-gene identifier including allele.

## Value

Character string. Lowercase nucleotide sequence starting at the
conserved C codon. Returns `""` if the gene is not found.
