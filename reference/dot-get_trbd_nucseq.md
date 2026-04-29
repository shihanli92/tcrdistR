# Get TRBD nucleotide sequences for an organism

Builds and caches a named list mapping 1-based integer D-gene IDs to
their nucleotide sequences.

## Usage

``` r
.get_trbd_nucseq(organism)
```

## Arguments

- organism:

  Character string. Organism identifier.

## Value

A named list with integer keys (as character names) and lowercase
nucleotide sequence values.
