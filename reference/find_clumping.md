# Find TCR clumping in TCR space

Identifies neighborhoods in TCR space containing more TCRs than expected
by chance under a null model of independent VDJ rearrangement. For each
TCR, counts how many other TCRs fall within a set of fixed TCRdist
radii, and compares the observed count to the Poisson expectation
derived from background distributions.

## Usage

``` r
find_clumping(
  tcr_df,
  organism,
  radii = c(24L, 48L, 72L, 96L),
  num_random_samples = 50000L,
  pvalue_threshold = 1,
  verbose = TRUE,
  clusters_gex = NULL,
  bg_tcrs = NULL,
  preserve_vj_pairings = FALSE
)
```

## Arguments

- tcr_df:

  A `data.frame` with columns `va`, `ja`, `cdr3a`, `cdr3a_nucseq`, `vb`,
  `jb`, `cdr3b`, `cdr3b_nucseq`.

- organism:

  Character string. Organism key (e.g. `"human"`, `"mouse"`).

- radii:

  Integer vector. TCRdist radii to test. Default
  `c(24L, 48L, 72L, 96L)`.

- num_random_samples:

  Integer. Number of random background chains per chain type. Default
  `50000L`.

- pvalue_threshold:

  Numeric. Maximum adjusted p-value to include in results. Default `1.0`
  (include all).

- verbose:

  Logical. Print progress messages. Default `TRUE`.

- clusters_gex:

  Integer vector of length `nrow(tcr_df)`, or `NULL`. If provided, also
  tests for TCR clumps within each GEX cluster. Default `NULL`.

- bg_tcrs:

  Optional data.frame. If provided, used for background generation
  instead of `tcr_df`. Default `NULL`.

- preserve_vj_pairings:

  Logical. Preserve V-J pairings in background resampling. Default
  `FALSE`.

## Value

A list with four elements:

- `results_df`:

  A data.frame sorted by `pvalue_adj` with columns: `clump_type`,
  `clone_index` (0-based), `nbr_radius`, `pvalue_adj`, `num_nbrs`,
  `expected_num_nbrs`, `raw_count`, `va`, `ja`, `cdr3a`, `vb`, `jb`,
  `cdr3b`, `clumping_group`, `clonotype_fdr_value`.

- `is_clumped`:

  Logical vector of length `nrow(tcr_df)`.

- `clusters`:

  Integer vector of length `nrow(tcr_df)`. 0 = not clumped, positive =
  cluster ID.

- `all_raw_pvalues`:

  Numeric matrix (`nrow(tcr_df)` x `length(radii)`).

## Details

The pipeline:

1.  Estimate per-TCR background frequency distributions via shuffled
    chain resampling.

2.  Assign alpha/beta chain groups for same-chain masking.

3.  Find all neighbors within `max(radii)` using
    [`tcrdist_radius_neighbors`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_radius_neighbors.md).

4.  Run Poisson tests at each radius (C++ via `rcpp_poisson_test_loop`).

5.  Perform single-linkage clustering of significant clumps.

## See also

[`setup_tcr_groups`](https://shihanli92.github.io/tcrdistR/reference/setup_tcr_groups.md),
[`tcrdist_radius_neighbors`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_radius_neighbors.md),
[`find_meta_clonotypes`](https://shihanli92.github.io/tcrdistR/reference/find_meta_clonotypes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- find_clumping(tcr_df, "human")
result$results_df
sum(result$is_clumped)
} # }
```
