
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tcrdistR

C++-accelerated TCR distance calculations for T-cell receptor repertoire
analysis. Computes pairwise TCRdist distances incorporating V-region and
CDR3 sequence comparisons using BLOSUM62-derived substitution matrices.
Targets feature parity with Python
[tcrdist3](https://github.com/kmayerb/tcrdist3).

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("shihanli92/tcrdistR")
```

## Quick Start

``` r
library(tcrdistR)

# Create a TCR data.frame (paired alpha/beta chains)
tcrs <- data.frame(
  va    = c("TRAV1-1*01", "TRAV1-2*01", "TRAV10*01",
            "TRAV1-1*01", "TRAV1-2*01"),
  cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF",
            "CAVRDSSYKLIF", "CAVKDSYKLIF"),
  vb    = c("TRBV5-1*01", "TRBV6-1*01", "TRBV7-2*01",
            "TRBV5-1*01", "TRBV6-1*01"),
  cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY",
            "CASSIRSSYEQF", "CASSIKSSYEQF"),
  stringsAsFactors = FALSE
)

# Compute the full pairwise distance matrix
dist_mat <- tcrdist_matrix(tcrs, "human")

# Find K-nearest neighbors
knn <- tcrdist_knn(tcrs, "human", K = 3)

# Sparse distance matrix (only store distances <= threshold)
sparse_mat <- tcrdist_sparse(tcrs, "human", threshold = 100)
```

## Features

**Distance computation** – Dense, sparse, and rectangular pairwise
TCRdist matrices with C++ backend: `tcrdist_matrix()`,
`tcrdist_sparse()`, `tcrdist_rect()`, `weighted_cdr3_distance()`,
`hamming_distance()`.

**Neighbor search** – K-nearest neighbors and radius-based queries:
`tcrdist_knn()`, `tcrdist_radius_neighbors()`, `knn_from_matrix()`,
`knn_from_pca()`.

**Clumping detection** – Background-corrected TCR clumping with Poisson
significance testing: `find_clumping()`, `setup_tcr_groups()`.

**Database matching** – Match query TCRs against literature epitope
databases with adjusted p-values: `find_significant_tcrdist_matches()`,
`match_tcrs_to_db()`.

**Dimensionality reduction** – Kernel PCA on TCRdist matrices (validated
against scipy): `compute_tcrdist_kernel_pca()`.

**CD8 scoring** – Logistic regression phenotype scoring:
`make_cd8_score_table_column()`.

**Diversity metrics** – Generalized Simpson’s entropy, fuzzy TCR-aware
diversity, richness, clonality: `tcr_diversity()`,
`tcr_fuzzy_diversity()`, `tcr_richness()`, `tcr_clonality()`.

**Clustering and neighborhood tests** – Hierarchical clustering with
Fisher’s exact / chi-squared neighborhood tests: `tcrdist_hclust()`,
`cluster_tcrs()`, `neighborhood_test()`.

**Meta-clonotypes** – Quasi-public TCR motif discovery across subjects:
`find_meta_clonotypes()`, `summarize_meta_clonotype()`.

**Distance-based joins** – Fuzzy join of TCR datasets by distance
threshold: `tcrdist_join()`.

**Visualization** – Distance heatmaps, dendrograms, CDR3 logos, gene
usage plots, PCA scatter: `plot_tcrdist_heatmap()`,
`plot_tcrdist_dendrogram()`, `plot_cdr3_logo()`, `plot_tcr_scatter()`,
and more.

**I/O** – Read TCR data from AIRR, Adaptive ImmunoSeq, 10x Genomics, and
generic CSV/TSV formats: `read_airr()`, `read_adaptive()`, `read_10x()`,
`read_tcr_table()`.

## Reading TCR Data

``` r
# Auto-detect format from common column naming conventions
tcrs <- read_tcr_table("my_tcr_data.csv")

# Format-specific readers
tcrs <- read_10x("filtered_contig_annotations.csv")
tcrs <- read_airr("airr_rearrangements.tsv")
tcrs <- read_adaptive("immunoseq_export.tsv")
```

## Documentation

- `vignette("tcrdistR-getting-started")` – Tutorial introduction
- `vignette("tcrdistR-advanced")` – Clumping, meta-clonotypes,
  clustering, DB matching
- `vignette("tcrdistR-visualization")` – Visualization guide

## License

MIT
