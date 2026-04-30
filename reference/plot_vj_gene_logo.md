# Plot a V/J gene usage logo

Renders gene names as stacked text glyphs with heights proportional to
their frequency — a "sequence logo" style for gene usage. Each gene name
is rendered as a coloured raster image and vertically stacked so that
more frequent genes occupy more vertical space.

## Usage

``` r
plot_vj_gene_logo(
  genes,
  organism,
  gene_type = c("V", "J"),
  chain = c("alpha", "beta"),
  max_genes = 10L
)
```

## Arguments

- genes:

  Character vector. V-gene or J-gene allele names (e.g.,
  `"TRAV1-2*01"`).

- organism:

  Character string. Organism identifier (e.g., `"human"`, `"mouse"`).

- gene_type:

  Character string. `"V"` or `"J"`.

- chain:

  Character string. `"alpha"` or `"beta"`.

- max_genes:

  Integer. Maximum number of genes to display. Default `10L`.

## Value

A `ggplot` object.

## See also

[`plot_gene_usage`](https://shihanli92.github.io/tcrdistR/reference/plot_gene_usage.md),
[`plot_tcr_logo_panel`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_logo_panel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(dash)
pa <- dash[dash$epitope == "PA", ]
plot_vj_gene_logo(pa$vb, organism = "mouse", gene_type = "V",
                  chain = "beta")
plot_vj_gene_logo(pa$ja, organism = "mouse", gene_type = "J",
                  chain = "alpha")
} # }
```
