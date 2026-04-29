# Create a TCRrep object

Constructs a `TCRrep` S4 object from a clonotype data frame and optional
parameters. This is the recommended way to create a `TCRrep` instance
(following Bioconductor convention, rather than calling
`new("TCRrep", ...)` directly).

## Usage

``` r
TCRrep(
  clone_df,
  organism = "human",
  chains = "AB",
  deduplicate = TRUE,
  metric = "tcrdist",
  compute_distances = FALSE,
  weight_cdr3 = WEIGHT_CDR3_REGION,
  gap_penalty_cdr3 = GAP_PENALTY_CDR3_REGION,
  weight_v_region = WEIGHT_V_REGION,
  gap_penalty_v_region = GAP_PENALTY_V_REGION
)
```

## Arguments

- clone_df:

  A `data.frame` of clonotypes. Required columns depend on the `chains`
  argument:

  `"AB"`

  :   Requires `va`, `cdr3a`, `vb`, `cdr3b`.

  `"A"`

  :   Requires `va`, `cdr3a`.

  `"B"`

  :   Requires `vb`, `cdr3b`.

  `"GD"`

  :   Requires `va`, `cdr3a`, `vb`, `cdr3b`.

- organism:

  Character string. Organism key recognised by
  [`load_gene_database`](https://shihanli92.github.io/tcrdistR/reference/load_gene_database.md),
  e.g. `"human"` or `"mouse"`.

- chains:

  Character string. One of `"AB"` (default), `"A"`, `"B"`, or `"GD"`.

- deduplicate:

  Controls clone deduplication (matching tcrdist3 behavior).

  `TRUE` (default)

  :   Deduplicate using chain columns (`va`, `cdr3a`, `vb`, `cdr3b`)
      plus `subject` if present. Within-subject duplicates are merged
      and `count` values summed.

  `FALSE`

  :   No deduplication; `clone_df` is stored as-is.

  Character vector

  :   Custom grouping columns. Only rows identical across all specified
      columns are merged. Example: `c("va", "cdr3a", "vb", "cdr3b")` to
      ignore subject.

- metric:

  Character string. Distance metric to use. One of `"tcrdist"`
  (default), `"hamming"`, `"levenshtein"`, or `"nw"`.

- compute_distances:

  Logical. If `TRUE` and `nrow(clone_df) > 0`, compute the pairwise
  distance matrix immediately and store it in the `paired_dist` slot.
  Defaults to `FALSE`.

- weight_cdr3:

  Integer. Weight applied to CDR3 distances. Defaults to
  `WEIGHT_CDR3_REGION` (3L).

- gap_penalty_cdr3:

  Integer. Gap penalty for CDR3 alignments. Defaults to
  `GAP_PENALTY_CDR3_REGION` (12L).

- weight_v_region:

  Integer. Weight applied to V-region distances. Defaults to
  `WEIGHT_V_REGION` (1L).

- gap_penalty_v_region:

  Integer. Gap penalty for V-region alignments. Defaults to
  `GAP_PENALTY_V_REGION` (4L).

## Value

A valid `TCRrep` S4 object.

## Details

Factor columns (`va`, `cdr3a`, `vb`, `cdr3b`) are automatically coerced
to character. The organism is validated against the bundled gene
database on construction.

## See also

[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md),
[`read_tcr_table`](https://shihanli92.github.io/tcrdistR/reference/read_tcr_table.md)

## Examples

``` r
# \donttest{
tcrs <- data.frame(
    va    = c("TRAV1-1*01", "TRAV1-1*01"),
    cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
    vb    = c("TRBV19*01", "TRBV19*01"),
    cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
    stringsAsFactors = FALSE
)

# Basic construction (deduplicates by default)
obj <- TCRrep(tcrs, organism = "human")

# With distance computation
obj <- TCRrep(tcrs, organism = "human", compute_distances = TRUE)
dim(obj@paired_dist)  # 2 x 2
#> [1] 2 2

# Custom dedup columns (ignore subject, collapse across individuals)
obj <- TCRrep(tcrs, organism = "human",
              deduplicate = c("va", "cdr3a", "vb", "cdr3b"))

# No deduplication
obj <- TCRrep(tcrs, organism = "human", deduplicate = FALSE)
# }
```
