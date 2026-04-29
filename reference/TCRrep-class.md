# TCRrep S4 class

Represents a T-cell receptor repertoire with associated metadata,
distance matrices, and analysis results. This is the central data object
of the tcrdistR package. Create instances using the
[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md)
constructor rather than calling `new("TCRrep", ...)` directly.

## Slots

- `clone_df`:

  A `data.frame` of clonotypes. Required columns depend on the `chains`
  value: alpha-beta (`"AB"`) requires `va`, `cdr3a`, `vb`, `cdr3b`;
  alpha only (`"A"`) requires `va`, `cdr3a`; beta only (`"B"`) requires
  `vb`, `cdr3b`; gamma-delta (`"GD"`) requires `va`, `cdr3a`, `vb`,
  `cdr3b`.

- `organism`:

  Character string specifying the organism. Must be one of the organisms
  available in the bundled gene database (e.g. `"human"`, `"mouse"`).

- `chains`:

  Character string specifying the chain combination. One of `"AB"`,
  `"A"`, `"B"`, or `"GD"`.

- `metric`:

  Character string specifying the distance metric. One of `"tcrdist"`,
  `"hamming"`, `"levenshtein"`, or `"nw"`.

- `weights`:

  A named list with elements `cdr3` and `v_region` specifying the
  integer weights applied to each region's distance contribution.

- `gap_penalties`:

  A named list with elements `cdr3` and `v_region` specifying the
  integer gap penalties used in sequence alignment.

- `dist_a`:

  Distance matrix for alpha-chain only comparisons. Either `NULL` or a
  numeric matrix.

- `dist_b`:

  Distance matrix for beta-chain only comparisons. Either `NULL` or a
  numeric matrix.

- `paired_dist`:

  Paired (alpha + beta combined) distance matrix. Either `NULL` or a
  numeric matrix. Populated when `compute_distances = TRUE` is passed to
  [`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md).

- `knn_indices`:

  K-nearest-neighbour index matrix. Either `NULL` or an integer matrix
  of dimensions N x K.

- `knn_distances`:

  K-nearest-neighbour distance matrix. Either `NULL` or a numeric matrix
  of dimensions N x K.

- `rep_dists`:

  A named list of representative-clonotype distance objects, populated
  by downstream analyses.

- `meta_clonotypes`:

  Meta-clonotype table. Either `NULL` or a `data.frame`.

- `motifs`:

  A named list of motif objects, populated by downstream analyses.

## See also

[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md)
for the recommended constructor.
