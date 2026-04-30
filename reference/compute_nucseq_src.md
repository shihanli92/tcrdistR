# Compute nucleotide source annotations for CDR3 sequences

For each TCR, calls
[`.analyze_junction()`](https://shihanli92.github.io/tcrdistR/reference/dot-analyze_junction.md)
to determine the V/N/D/J origin of each nucleotide in the CDR3 region.
The result can be passed to
[`plot_cdr3_logo`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md)
via its `nucseq_src` parameter to display junction bars showing the
rearrangement structure.

## Usage

``` r
compute_nucseq_src(tcrs, organism, chain = c("alpha", "beta"))
```

## Arguments

- tcrs:

  Data.frame with columns `va`, `ja`, `cdr3a`, `cdr3a_nucseq` (for alpha
  chain) or `vb`, `jb`, `cdr3b`, `cdr3b_nucseq` (for beta chain).

- organism:

  Character string. Organism identifier (e.g., `"human"`, `"mouse"`).

- chain:

  Character string. `"alpha"` or `"beta"`.

## Value

A list of character vectors (one per TCR). Each vector has length
`nchar(cdr3_nucseq)` with elements from `c("V", "N", "J")` for alpha
chain, or `c("V", "N1", "D", "N2", "J")` for beta chain. Returns `NULL`
for TCRs where junction analysis fails.

## See also

[`plot_cdr3_logo`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md),
[`plot_junction_bars`](https://shihanli92.github.io/tcrdistR/reference/plot_junction_bars.md),
[`plot_tcr_logo_panel`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_logo_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(dash)
pa <- dash[dash$epitope == "PA", ][1:20, ]
src <- compute_nucseq_src(pa, organism = "mouse", chain = "beta")
plot_cdr3_logo(pa$cdr3b, chain = "beta", nucseq_src = src,
               show_junction_bars = TRUE)
} # }
```
