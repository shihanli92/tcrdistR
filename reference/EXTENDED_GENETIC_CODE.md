# Extended genetic code including degenerate IUPAC codons

A named character vector extending `GENETIC_CODE` with all degenerate
IUPAC codons formed from the 15 symbols in `NUCLEOTIDE_CLASSES`. For
each degenerate codon, all possible standard expansions are looked up in
`GENETIC_CODE`:

- If all expansions map to the same amino acid, that amino acid is
  assigned.

- If expansions map to different amino acids, `"X"` is assigned.

Standard 64 codons retain their original mappings. Total size is
`15^3 = 3375` entries (all possible combinations of the 15 IUPAC
symbols).

## Usage

``` r
EXTENDED_GENETIC_CODE
```

## Format

A named character vector of length 3375.

## Details

Built once at package load time via
[`local()`](https://rdrr.io/r/base/eval.html).
