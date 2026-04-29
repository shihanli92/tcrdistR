# Load and cache the TCR/BCR gene database

Parses the bundled gene database TSV on first call, computes sequence
representatives (exact and mm1) for V and J genes, and caches the result
in the package-private environment `.tcrdistR_env$all_genes`. Subsequent
calls return the cached data without re-parsing.

## Usage

``` r
load_gene_database(organism = NULL)
```

## Arguments

- organism:

  Character string or `NULL`. If a non-`NULL` string is supplied, only
  the gene entries for that organism are returned. If `NULL` (default),
  the complete nested list for all organisms is returned.

## Value

When `organism` is `NULL`: a named list keyed by organism name, where
each element is itself a named list of gene entries keyed by gene ID.

When `organism` is a character string: a named list of gene entries for
that organism, keyed by gene ID.

Each gene entry is a named list with fields:

- id:

  Character. Gene identifier including allele, e.g. `"TRAV1*01"`.

- organism:

  Character. Organism name.

- chain:

  Character. Either `"A"` (alpha / gamma) or `"B"` (beta / delta).

- region:

  Character. Gene segment: `"V"`, `"D"`, or `"J"`.

- nucseq:

  Character. Nucleotide sequence.

- alseq:

  Character. Aligned protein sequence (gaps represented as `"."`).

- cdrs:

  Character vector. CDR subsequences extracted from `alseq`.

- cdr_columns:

  List of 2-element integer vectors. Start and end positions (1-indexed,
  inclusive) of each CDR in `alseq`.

- nucseq_offset:

  Integer. 0-based reading frame offset (frame - 1).

- protseq:

  Character. Protein sequence without gap characters.

- rep:

  Character. Representative gene ID from exact loopseq neighbours (min
  by ID).

- mm1_rep:

  Character. Representative gene ID from transitive mm1 loopseq
  neighbours.

- count_rep:

  Character. Gene-level name (allele stripped); used for clone counting.

## Details

The database file is `inst/extdata/combo_xcr_2023-12-30.tsv` and
contains 2836 gene entries across organisms including `"human"`,
`"mouse"`, `"human_ig"`, `"mouse_ig"`, `"human_gd"`, `"mouse_gd"`, and
`"rhesus"`.

## See also

[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md),
[`trim_allele_to_gene`](https://shihanli92.github.io/tcrdistR/reference/trim_allele_to_gene.md)

## Examples

``` r
# \donttest{
# Load all organisms
all_g <- load_gene_database()
names(all_g)  # "human", "mouse", ...
#> [1] "mouse"     "human"     "mouse_gd"  "human_gd"  "rhesus"    "rhesus_gd"
#> [7] "mouse_ig"  "human_ig" 

# Load one organism
human_genes <- load_gene_database("human")
human_genes[["TRAV1-1*01"]]$protseq
#> [1] "GQSLEQPSEVTAVEGAIVQINCTYQTSGFYGLSWYQQHDGGAPTFLSYNALDGLEETGRFSSFLSRSDSYGYLLLQELQMKDSASYFCAVR"
# }
```
