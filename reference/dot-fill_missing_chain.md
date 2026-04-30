# Fill missing chain columns for single-chain input

Detects whether both chains are present and, if not, fills the missing
chain with dummy values so that downstream code can proceed unchanged.
The dummy V-gene is the first V-gene for that chain in the gene
database; the dummy CDR3 is a minimal valid sequence. When combined with
the `components` parameter (auto-set to `"alpha"` or `"beta"`), the
dummy chain contributes exactly zero to the distance.

## Usage

``` r
.fill_missing_chain(tcrs, organism)
```

## Arguments

- tcrs:

  A `data.frame`.

- organism:

  Character string passed to `load_gene_database`.

## Value

A list with elements `tcrs` (possibly augmented) and `chain` (`"AB"`,
`"A"`, or `"B"`).
