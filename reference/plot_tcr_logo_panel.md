# Plot a composite TCR rearrangement logo panel

Arranges V-gene logos, CDR3 sequence logos (with optional junction
bars), and J-gene logos for both alpha and beta chains in a single
composite panel. This provides a comprehensive view of TCR rearrangement
structure.

## Usage

``` r
plot_tcr_logo_panel(tcrs, organism, show_junction_bars = TRUE, title = NULL)
```

## Arguments

- tcrs:

  Data.frame with at least columns `va`, `ja`, `cdr3a`, `vb`, `jb`,
  `cdr3b`. For junction bars, also requires `cdr3a_nucseq` and
  `cdr3b_nucseq`.

- organism:

  Character string. Organism identifier (e.g., `"human"`, `"mouse"`).

- show_junction_bars:

  Logical. If `TRUE` (default) and nucleotide sequence columns are
  present, display V/N/D/J junction bars below each CDR3 logo.

- title:

  Optional character string. Overall panel title.

## Value

A `patchwork` object.

## See also

[`plot_vj_gene_logo`](https://shihanli92.github.io/tcrdistR/reference/plot_vj_gene_logo.md),
[`plot_cdr3_logo`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md),
[`plot_junction_bars`](https://shihanli92.github.io/tcrdistR/reference/plot_junction_bars.md),
[`compute_nucseq_src`](https://shihanli92.github.io/tcrdistR/reference/compute_nucseq_src.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(dash)
pa <- dash[dash$epitope == "PA", ][1:30, ]
plot_tcr_logo_panel(pa, organism = "mouse")
plot_tcr_logo_panel(pa, organism = "mouse", show_junction_bars = FALSE)
} # }
```
