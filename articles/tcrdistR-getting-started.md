# Getting Started with tcrdistR

## Introduction

**tcrdistR** computes pairwise distances between T-cell receptors (TCRs)
using the TCRdist metric. Each TCR is represented by its V-gene and CDR3
amino acid sequence for both alpha and beta chains. The distance
incorporates V-region similarity (via BLOSUM62-derived substitution
matrices) and CDR3 sequence alignment with variable gap positioning
matching tcrdist3.

All distance computations are implemented in C++ for high performance.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("shihanli92/tcrdistR")
```

## The DASH Dataset

tcrdistR ships with the **DASH** dataset (Dash et al., 2017) — 1924
paired alpha-beta mouse TCRs responding to 7 viral epitopes, collected
from 78 subjects. We use it throughout these vignettes.

``` r
library(tcrdistR)
data(dash)

dim(dash)
#> [1] 1924   12
colnames(dash)
#>  [1] "subject"      "epitope"      "count"        "va"           "ja"          
#>  [6] "cdr3a"        "cdr3a_nucseq" "vb"           "jb"           "cdr3b"       
#> [11] "cdr3b_nucseq" "clone_id"
table(dash$epitope)
#> 
#>   F2 m139  M38  M45   NP   PA  PB1 
#>  117   87  158  291  305  324  642
```

The required columns for distance computation are `va`, `cdr3a`, `vb`,
and `cdr3b`. Additional columns like `epitope` and `subject` are carried
along for downstream analysis.

## Reading TCR Data from Files

tcrdistR can read TCR data from multiple common formats. All readers
produce a data.frame with standardized column names (`va`, `cdr3a`,
`vb`, `cdr3b`).

``` r
# Auto-detect format from column naming conventions
tcrs <- read_tcr_table("my_data.tsv")

# Format-specific readers
tcrs <- read_10x("filtered_contig_annotations.csv")
tcrs <- read_airr("airr_rearrangements.tsv")
tcrs <- read_adaptive("immunoseq_export.tsv")
```

## Computing the Distance Matrix

The core function
[`tcrdist_matrix()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)
computes a dense N x N pairwise distance matrix. We use a subset of 50
PA-specific TCRs for speed:

``` r
pa <- dash[dash$epitope == "PA", ]
pa_sub <- pa[1:50, ]

dist_mat <- tcrdist_matrix(pa_sub, organism = "mouse")
dim(dist_mat)
#> [1] 50 50
dist_mat[1:5, 1:5]
#>     1   2   3   4   5
#> 1   0 317 341 326 281
#> 2 317   0 177 243 204
#> 3 341 177   0 288 216
#> 4 326 243 288   0 138
#> 5 281 204 216 138   0
```

Each entry is an integer TCRdist value. Identical TCRs have distance 0;
typical distances range from 0 to ~400.

## Sparse Distances

For large datasets, storing the full N x N matrix is impractical. Use
[`tcrdist_sparse()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_sparse.md)
to only store distances below a threshold:

``` r
sparse_mat <- tcrdist_sparse(pa_sub, organism = "mouse", threshold = 100)
class(sparse_mat)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
```

The result is a `Matrix::dgCMatrix`. Non-stored entries represent
distances above the threshold.

## Rectangular Distances

Compute distances between two different sets of TCRs (query vs
reference):

``` r
query <- pa_sub[1:5, ]
reference <- pa_sub[6:50, ]
rect_mat <- tcrdist_rect(query, reference, organism = "mouse")
dim(rect_mat)
#> [1]  5 45
```

## K-Nearest Neighbors

Find the K closest TCRs for each input. This computes distances and
extracts neighbors in a single pass without materializing the full N x N
matrix:

``` r
knn <- tcrdist_knn(pa_sub, organism = "mouse", K = 5L)
dim(knn$knn_indices)    # N x K
#> [1] 50  5
dim(knn$knn_distances)  # N x K
#> [1] 50  5

# Nearest neighbor of TCR 1
knn$knn_indices[1, ]
#> [1] 49 43  9 14 15
knn$knn_distances[1, ]
#> [1] 211 227 235 245 248
```

## Radius-Based Neighbors

Find all neighbors within a distance threshold:

``` r
neighbors <- tcrdist_radius_neighbors(pa_sub, organism = "mouse", radius = 50)
length(neighbors)  # one entry per TCR
#> [1] 50

# First TCR's neighbors
neighbors[[1]]$indices
#> integer(0)
neighbors[[1]]$distances
#> numeric(0)
```

## Per-Component Distances

All distance functions accept a `components` parameter to compute
distances from specific TCR regions. This is useful for comparing
CDR3-only distances, V-region-only contributions, or single-chain
distances.

``` r
# CDR3 only (both chains, no V-region)
d_cdr3 <- tcrdist_matrix(pa_sub, "mouse", components = "cdr3")

# V-region only (CDR1 + CDR2 + CDR2.5)
d_v <- tcrdist_matrix(pa_sub, "mouse", components = "v_region")

# Components are additive: v_region + cdr3 = all
d_all <- tcrdist_matrix(pa_sub, "mouse")
stopifnot(all(abs((d_v + d_cdr3) - d_all) < 1e-10))

# Single chain
d_alpha <- tcrdist_matrix(pa_sub, "mouse", components = "alpha")
d_beta  <- tcrdist_matrix(pa_sub, "mouse", components = "beta")
stopifnot(all(abs((d_alpha + d_beta) - d_all) < 1e-10))

# Custom combination
d_custom <- tcrdist_matrix(pa_sub, "mouse", components = c("cdr3a", "vb"))
```

Available presets: `"all"` (default), `"cdr3"`, `"v_region"`, `"alpha"`,
`"beta"`. Individual terms: `"va"`, `"cdr3a"`, `"vb"`, `"cdr3b"`.

## Single-Chain Mode

Many datasets contain only beta-chain data (e.g., from bulk sequencing).
All distance functions automatically detect single-chain input and
compute distances for the available chain only:

``` r
# Beta-only data (no alpha columns)
beta_only <- pa_sub[, c("vb", "cdr3b")]
d_beta_only <- tcrdist_matrix(beta_only, "mouse")
dim(d_beta_only)
#> [1] 50 50

# Result matches components="beta" on paired data
d_beta_paired <- tcrdist_matrix(pa_sub, "mouse", components = "beta")
stopifnot(all(abs(d_beta_only - d_beta_paired) < 1e-10))

# Works with all wrappers: sparse, rect, knn, radius_neighbors
knn_beta <- tcrdist_knn(beta_only, "mouse", K = 3L)
```

## Standardizing TCR Data

Use
[`as_tcr_df()`](https://shihanli92.github.io/tcrdistR/reference/as_tcr_df.md)
to convert data.frames from common tools (scRepertoire, scirpy,
dandelion, tcrdist3) to tcrdistR’s canonical column names:

``` r
# Auto-detect format
tcrs <- as_tcr_df(my_data)

# Explicit format
tcrs <- as_tcr_df(my_data, format = "screpertoire")

# Custom column mapping
tcrs <- as_tcr_df(my_data, col_map = c(va = "V_alpha", cdr3a = "CDR3_alpha",
                                        vb = "V_beta",  cdr3b = "CDR3_beta"))
```

## Single-Pair Distances

For comparing individual CDR3 sequences:

``` r
# Weighted CDR3 distance (uses BLOSUM62-derived BSD4 matrix)
weighted_cdr3_distance("CAVSLDSNYQLIW", "CALGDRATGGNNKLTF")
#> [1] 102

# Simple Hamming distance (counts mismatches at each position)
hamming_distance("CASSI", "CASSK")
#> [1] 1
```

## Diversity Metrics

Repertoire diversity functions work on clone count vectors:

``` r
# Clone counts for PA-specific TCRs
pa_counts <- pa$count

# Generalized diversity (order 2 = inverse Simpson)
div <- tcr_diversity(pa_counts, order = 2)
div$entropy
#> [1] 0.9958049
div$effective_number
#> [1] 238.3727

# Species richness
tcr_richness(pa_counts)
#> [1] 324

# Clonality (1 - normalized Shannon entropy)
tcr_clonality(pa_counts)
#> [1] 0.05175504
```

## The TCRrep Object

For more structured workflows, wrap your data in a `TCRrep` S4 object.
By default,
[`TCRrep()`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md)
deduplicates identical clones (matching tcrdist3 behavior), grouping by
chain columns plus `subject` if present, and summing clone counts:

``` r
rep <- TCRrep(pa_sub, organism = "mouse")
#> deduplicate: 50 -> 49 clones
nrow(pa_sub)        # before dedup
#> [1] 50
nrow(rep@clone_df)  # after dedup (within-subject duplicates merged)
#> [1] 49
rep
#> TCRrep object: 49 clonotypes
#>   organism: mouse
#>   chains: AB
#>   metric: tcrdist
#>   distances: not computed
```

You can control deduplication with the `deduplicate` argument:

``` r
# Chain columns only (collapse across subjects)
rep_strict <- TCRrep(pa_sub, organism = "mouse", deduplicate = character(0))
#> deduplicate: 50 -> 47 clones
nrow(rep_strict@clone_df)
#> [1] 47

# No deduplication
rep_raw <- TCRrep(pa_sub, organism = "mouse", deduplicate = FALSE)
nrow(rep_raw@clone_df)
#> [1] 50
```

## Next Steps

- [`vignette("tcrdistR-advanced")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-advanced.md)
  — Clustering, neighborhood tests, meta-clonotypes, database matching,
  kernel PCA, fuzzy diversity
- [`vignette("tcrdistR-visualization")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-visualization.md)
  — Heatmaps, dendrograms, CDR3 logos, scatter plots
