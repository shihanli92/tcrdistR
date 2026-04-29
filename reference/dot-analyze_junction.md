# Analyze a TCR junction region

Aligns the CDR3 nucleotide sequence against germline V and J sequences
to identify trimming, N-insertion lengths, and D-gene usage (beta only).

## Usage

``` r
.analyze_junction(
  organism,
  v_gene,
  j_gene,
  cdr3_protseq,
  cdr3_nucseq,
  force_d_id = 0L,
  mismatch_score = .DEFAULT_MISMATCH_SCORE_JUNCTION_ANALYSIS
)
```

## Arguments

- organism:

  Character string. Organism identifier.

- v_gene:

  Character string. V-gene identifier including allele.

- j_gene:

  Character string. J-gene identifier including allele.

- cdr3_protseq:

  Character string. CDR3 amino acid sequence.

- cdr3_nucseq:

  Character string. CDR3 nucleotide sequence (lowercase).

- force_d_id:

  Integer. If non-zero, only consider this D-gene ID.

- mismatch_score:

  Integer. Score for a mismatch.

## Value

A named list with elements: `new_nucseq`, `cdr3_protseq_masked`,
`cdr3_protseq_new_nuc_countstring`, `cdr3_nucseq_src`, `trims`,
`inserts`.
