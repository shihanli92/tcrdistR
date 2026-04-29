# Translate a nucleotide sequence to a protein sequence

Translates a nucleotide sequence into a single-letter amino acid string
using `EXTENDED_GENETIC_CODE`. Supports all IUPAC degenerate nucleotide
codes. Codons containing `"#"` are translated as `"#"` (gap indicator
used in some CoNGA input files). Codons not found in the extended code
are translated as `"X"`.

## Usage

``` r
get_translation(seq, frame = "+1")
```

## Arguments

- seq:

  Character string. The nucleotide sequence to translate. May contain
  standard or IUPAC-degenerate bases and may be any case (converted to
  lowercase internally after offset trimming).

- frame:

  Character string of length 1. Reading frame, e.g. `"+1"`, `"+2"`,
  `"+3"` for forward strand, or `"-1"`, `"-2"`, `"-3"` for
  reverse-complement strand. The sign determines strand and the absolute
  integer value (1, 2, or 3) determines the 0-based offset into the
  sequence before translation begins.

## Value

A character string of amino acids. Length is
`floor((nchar(seq) - offset) / 3)`.

## Examples

``` r
# Translate from frame +1 (no offset)
get_translation("ATGAAATTT", "+1")  # "MKF"
#> [1] "MKF"

# Frame +2 skips the first nucleotide
get_translation("AATGAAATTT", "+2")  # "MKF"
#> [1] "MKF"
```
