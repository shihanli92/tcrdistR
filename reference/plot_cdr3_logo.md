# Plot a CDR3 sequence logo

Renders a sequence logo from CDR3 amino acid sequences using ggseqlogo.
Two methods are available: `"bits"` shows information content per
position, `"prob"` shows amino acid frequencies.

## Usage

``` r
plot_cdr3_logo(
  cdr3_seqs,
  chain = c("alpha", "beta"),
  trim = TRUE,
  method = c("prob", "bits"),
  gap_character = "-",
  title = NULL,
  nucseq_src = NULL,
  show_junction_bars = FALSE
)
```

## Arguments

- cdr3_seqs:

  Character vector. CDR3 amino acid sequences.

- chain:

  Character string. `"alpha"` or `"beta"`.

- trim:

  Logical. If `TRUE`, trim conserved C/F residues before alignment
  (first 3, last 2 characters). Default `TRUE`.

- method:

  Character string. `"prob"` for frequency-based logo, `"bits"` for
  information content. Default `"prob"`.

- gap_character:

  Character string. Gap character for alignment. Default `"-"`.

- title:

  Optional plot title. If `NULL`, defaults to the chain name with Greek
  letter.

- nucseq_src:

  List of character vectors with nucleotide source labels (one per
  CDR3). Required for junction bars.

- show_junction_bars:

  Logical. If `TRUE` and `nucseq_src` is provided, stack junction bars
  below the logo via patchwork.

## Value

A `ggplot` object (or `patchwork` object if junction bars are included).

## See also

[`plot_junction_bars`](https://shihanli92.github.io/tcrdistR/reference/plot_junction_bars.md),
[`plot_gene_usage`](https://shihanli92.github.io/tcrdistR/reference/plot_gene_usage.md)

## Examples

``` r
if (FALSE) { # \dontrun{
seqs <- c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF",
          "CAVRDSSYKLIF", "CAVKDSYKLIF")
plot_cdr3_logo(seqs, chain = "alpha", method = "bits")
plot_cdr3_logo(seqs, chain = "beta", method = "prob")
} # }
```
