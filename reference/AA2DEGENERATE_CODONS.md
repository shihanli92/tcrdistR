# Degenerate codon representations per amino acid

Named list mapping each single-letter amino acid code to a character
vector of IUPAC degenerate codon strings. Amino acids whose codons span
multiple two-letter prefixes have multiple entries (e.g., L = c("ctn",
"ttr")). Stop codons are excluded, matching the Python source.

## Usage

``` r
AA2DEGENERATE_CODONS
```

## Format

A named list of length 20.

## Details

The degeneracy mapping used for the third nucleotide position:

- a -\> a, c -\> c, g -\> g, t -\> t

- ct -\> y, ag -\> r

- act -\> h

- acgt -\> n
