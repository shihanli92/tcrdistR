# Plot junction bars showing V/N/D/J nucleotide composition

Renders stacked coloured bars showing the fraction of nucleotides
derived from V-gene, N-region, D-gene, or J-gene segments at each codon
position within a CDR3 alignment. Designed to sit below a CDR3 sequence
logo.

## Usage

``` r
plot_junction_bars(
  junction_pwm,
  chain = c("alpha", "beta"),
  gap_character = "-"
)
```

## Arguments

- junction_pwm:

  Numeric matrix. Rows are junction source labels, columns are
  nucleotide positions (3 per codon).

- chain:

  Character string. `"alpha"` or `"beta"`.

- gap_character:

  Character string. Gap character. Default `"-"`.

## Value

A `ggplot` object.

## See also

[`plot_cdr3_logo`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md)
