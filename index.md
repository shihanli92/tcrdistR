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
TCRdist matrices with C++ backend:
[`tcrdist_matrix()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md),
[`tcrdist_sparse()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_sparse.md),
[`tcrdist_rect()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_rect.md),
[`weighted_cdr3_distance()`](https://shihanli92.github.io/tcrdistR/reference/weighted_cdr3_distance.md),
[`hamming_distance()`](https://shihanli92.github.io/tcrdistR/reference/hamming_distance.md).

**Neighbor search** – K-nearest neighbors and radius-based queries:
[`tcrdist_knn()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_knn.md),
[`tcrdist_radius_neighbors()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_radius_neighbors.md),
[`knn_from_matrix()`](https://shihanli92.github.io/tcrdistR/reference/knn_from_matrix.md),
[`knn_from_pca()`](https://shihanli92.github.io/tcrdistR/reference/knn_from_pca.md).

**Clumping detection** – Background-corrected TCR clumping with Poisson
significance testing:
[`find_clumping()`](https://shihanli92.github.io/tcrdistR/reference/find_clumping.md),
[`setup_tcr_groups()`](https://shihanli92.github.io/tcrdistR/reference/setup_tcr_groups.md).

**Database matching** – Match query TCRs against literature epitope
databases with adjusted p-values:
[`find_significant_tcrdist_matches()`](https://shihanli92.github.io/tcrdistR/reference/find_significant_tcrdist_matches.md),
[`match_tcrs_to_db()`](https://shihanli92.github.io/tcrdistR/reference/match_tcrs_to_db.md).

**Dimensionality reduction** – Kernel PCA on TCRdist matrices (validated
against scipy):
[`compute_tcrdist_kernel_pca()`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_kernel_pca.md).

**CD8 scoring** – Logistic regression phenotype scoring:
[`make_cd8_score_table_column()`](https://shihanli92.github.io/tcrdistR/reference/make_cd8_score_table_column.md).

**Diversity metrics** – Generalized Simpson’s entropy, fuzzy TCR-aware
diversity, richness, clonality:
[`tcr_diversity()`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md),
[`tcr_fuzzy_diversity()`](https://shihanli92.github.io/tcrdistR/reference/tcr_fuzzy_diversity.md),
[`tcr_richness()`](https://shihanli92.github.io/tcrdistR/reference/tcr_richness.md),
[`tcr_clonality()`](https://shihanli92.github.io/tcrdistR/reference/tcr_clonality.md).

**Clustering and neighborhood tests** – Hierarchical clustering with
Fisher’s exact / chi-squared neighborhood tests:
[`tcrdist_hclust()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_hclust.md),
[`cluster_tcrs()`](https://shihanli92.github.io/tcrdistR/reference/cluster_tcrs.md),
[`neighborhood_test()`](https://shihanli92.github.io/tcrdistR/reference/neighborhood_test.md).

**Meta-clonotypes** – Quasi-public TCR motif discovery across subjects:
[`find_meta_clonotypes()`](https://shihanli92.github.io/tcrdistR/reference/find_meta_clonotypes.md),
[`summarize_meta_clonotype()`](https://shihanli92.github.io/tcrdistR/reference/summarize_meta_clonotype.md).

**Distance-based joins** – Fuzzy join of TCR datasets by distance
threshold:
[`tcrdist_join()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_join.md).

**Visualization** – Distance heatmaps, dendrograms, CDR3 logos, gene
usage plots, PCA scatter:
[`plot_tcrdist_heatmap()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_heatmap.md),
[`plot_tcrdist_dendrogram()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_dendrogram.md),
[`plot_cdr3_logo()`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md),
[`plot_tcr_scatter()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_scatter.md),
and more.

**I/O** – Read TCR data from AIRR, Adaptive ImmunoSeq, 10x Genomics, and
generic CSV/TSV formats:
[`read_airr()`](https://shihanli92.github.io/tcrdistR/reference/read_airr.md),
[`read_adaptive()`](https://shihanli92.github.io/tcrdistR/reference/read_adaptive.md),
[`read_10x()`](https://shihanli92.github.io/tcrdistR/reference/read_10x.md),
[`read_tcr_table()`](https://shihanli92.github.io/tcrdistR/reference/read_tcr_table.md).

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

- [`vignette("tcrdistR-getting-started")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-getting-started.md)
  – Tutorial introduction
- [`vignette("tcrdistR-advanced")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-advanced.md)
  – Clumping, meta-clonotypes, clustering, DB matching
- [`vignette("tcrdistR-visualization")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-visualization.md)
  – Visualization guide

## License

MIT
