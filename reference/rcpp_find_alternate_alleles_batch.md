# Find alternate alleles for a batch of TCRs (C++ implementation)

For each paired TCR (alpha+beta chains), tries alternate V/J alleles
(same gene, different allele number after `*`) and picks the allele
whose CDR3 nucleotide prefix match score improves by at least
`min_improvement` over the current allele. J alleles are compared using
reversed sequences.

## Usage

``` r
rcpp_find_alternate_alleles_batch(
  va_genes,
  ja_genes,
  cdr3a_nucseqs,
  vb_genes,
  jb_genes,
  cdr3b_nucseqs,
  all_gene_names,
  all_v_cdr3_nucseqs,
  all_j_cdr3_nucseqs,
  mismatch_score = -4L,
  min_improvement = 2L
)
```

## Arguments

- va_genes:

  Character vector. Alpha V gene names.

- ja_genes:

  Character vector. Alpha J gene names.

- cdr3a_nucseqs:

  Character vector. Alpha CDR3 nucleotide sequences.

- vb_genes:

  Character vector. Beta V gene names.

- jb_genes:

  Character vector. Beta J gene names.

- cdr3b_nucseqs:

  Character vector. Beta CDR3 nucleotide sequences.

- all_gene_names:

  Character vector. All gene names in the database.

- all_v_cdr3_nucseqs:

  Character vector. Pre-computed V CDR3 nucleotide sequences (same order
  as `all_gene_names`).

- all_j_cdr3_nucseqs:

  Character vector. Pre-computed J CDR3 nucleotide sequences (same order
  as `all_gene_names`).

- mismatch_score:

  Integer. Mismatch penalty for prefix scoring. Default -4.

- min_improvement:

  Integer. Minimum score improvement required to accept an alternate
  allele over the current one. Default 2.

## Value

A list with elements: `new_va`, `new_ja`, `new_vb`, `new_jb` (character
vectors), and `counts` (data.frame with columns `old_gene`, `new_gene`,
`count`).

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- rcpp_find_alternate_alleles_batch(
    va_genes, ja_genes, cdr3a_nucseqs,
    vb_genes, jb_genes, cdr3b_nucseqs,
    all_gene_names, all_v_cdr3_nucseqs, all_j_cdr3_nucseqs
  )
} # }
```
