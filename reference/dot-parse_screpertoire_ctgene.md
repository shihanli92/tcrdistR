# Parse scRepertoire CTgene column into V/J gene columns

Parses scRepertoire's concatenated gene format:
`TRAV.TRAJ.TRAC_TRBV.TRBD.TRBJ.TRBC` (alpha V.J.C underscore beta
V.D.J.C). Multi-chain cells (semicolon-separated) use the first chain
per locus.

## Usage

``` r
.parse_screpertoire_ctgene(ctgene_col)
```

## Arguments

- ctgene_col:

  Character vector. The `CTgene` column values.

## Value

A `data.frame` with columns `va`, `ja`, `vb`, `jb`.
