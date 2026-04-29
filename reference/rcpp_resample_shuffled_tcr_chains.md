# Resample shuffled TCR chains by junction breakpoint splicing (C++ implementation)

Randomly selects pairs of junctions, finds shared breakpoints, and
splices chimeric nucleotide sequences. Validates that the result has a
length divisible by 3 and contains no stop codons, then translates to
amino acid.

## Usage

``` r
rcpp_resample_shuffled_tcr_chains(
  junction_v_genes,
  junction_j_genes,
  junction_nucseqs,
  junction_breakpoints_pre_d,
  junction_breakpoints_post_d,
  chain,
  num_samples,
  max_attempts
)
```

## Arguments

- junction_v_genes:

  Character vector. V gene names for each junction.

- junction_j_genes:

  Character vector. J gene names for each junction.

- junction_nucseqs:

  Character vector. CDR3 nucleotide sequences.

- junction_breakpoints_pre_d:

  List of integer vectors. Pre-D breakpoints for each junction.

- junction_breakpoints_post_d:

  List of integer vectors. Post-D breakpoints for each junction.

- chain:

  Character string. `"A"` for alpha (pre-D breakpoints only) or `"B"`
  for beta (both pre-D and post-D).

- num_samples:

  Integer. Number of chimeric TCRs to generate.

- max_attempts:

  Integer. Maximum sampling attempts before returning partial results.

## Value

A data.frame with columns: `v_gene`, `j_gene`, `cdr3` (amino acid),
`cdr3_nucseq`.

## Details

Uses an inline standard genetic code table (no external codon map
needed).

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- rcpp_resample_shuffled_tcr_chains(
    v_genes, j_genes, nucseqs, bp_pre, bp_post, "B", 1000L, 100000L
  )
} # }
```
