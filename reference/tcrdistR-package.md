# tcrdistR: C++-Accelerated TCR Distance Calculations

Compute pairwise TCRdist distances for T-cell receptor (TCR) repertoire
analysis. All distance computations are implemented in C++ via Rcpp for
high performance, with a thin R layer for input validation, I/O, and
visualization. Targets feature parity with Python tcrdist3.

## Details

The package provides the following functional areas:

**Distance computation:**
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md),
[`tcrdist_sparse`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_sparse.md),
[`tcrdist_rect`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_rect.md),
[`weighted_cdr3_distance`](https://shihanli92.github.io/tcrdistR/reference/weighted_cdr3_distance.md),
[`bsd4_matrix`](https://shihanli92.github.io/tcrdistR/reference/bsd4_matrix.md),
[`hamming_distance`](https://shihanli92.github.io/tcrdistR/reference/hamming_distance.md),
[`hamming_matrix`](https://shihanli92.github.io/tcrdistR/reference/hamming_matrix.md)

**Neighbor search:**
[`tcrdist_knn`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_knn.md),
[`tcrdist_radius_neighbors`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_radius_neighbors.md),
[`knn_from_matrix`](https://shihanli92.github.io/tcrdistR/reference/knn_from_matrix.md),
[`knn_from_pca`](https://shihanli92.github.io/tcrdistR/reference/knn_from_pca.md)

**Clumping and background models:**
[`find_clumping`](https://shihanli92.github.io/tcrdistR/reference/find_clumping.md),
[`setup_tcr_groups`](https://shihanli92.github.io/tcrdistR/reference/setup_tcr_groups.md)

**Database matching:**
[`find_significant_tcrdist_matches`](https://shihanli92.github.io/tcrdistR/reference/find_significant_tcrdist_matches.md),
[`match_tcrs_to_db`](https://shihanli92.github.io/tcrdistR/reference/match_tcrs_to_db.md),
[`strict_single_chain_match_tcrs_to_db`](https://shihanli92.github.io/tcrdistR/reference/strict_single_chain_match_tcrs_to_db.md)

**Dimensionality reduction:**
[`compute_tcrdist_kernel_pca`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_kernel_pca.md)

**CD8 scoring:**
[`make_cd8_score_table_column`](https://shihanli92.github.io/tcrdistR/reference/make_cd8_score_table_column.md)

**Diversity metrics:**
[`tcr_diversity`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md),
[`tcr_fuzzy_diversity`](https://shihanli92.github.io/tcrdistR/reference/tcr_fuzzy_diversity.md),
[`tcr_richness`](https://shihanli92.github.io/tcrdistR/reference/tcr_richness.md),
[`tcr_clonality`](https://shihanli92.github.io/tcrdistR/reference/tcr_clonality.md)

**Clustering and neighborhood tests:**
[`tcrdist_hclust`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_hclust.md),
[`cluster_tcrs`](https://shihanli92.github.io/tcrdistR/reference/cluster_tcrs.md),
[`neighborhood_test`](https://shihanli92.github.io/tcrdistR/reference/neighborhood_test.md)

**Meta-clonotypes:**
[`find_meta_clonotypes`](https://shihanli92.github.io/tcrdistR/reference/find_meta_clonotypes.md),
[`summarize_meta_clonotype`](https://shihanli92.github.io/tcrdistR/reference/summarize_meta_clonotype.md)

**Joins:**
[`tcrdist_join`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_join.md)

**Visualization:**
[`plot_tcrdist_heatmap`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_heatmap.md),
[`plot_tcrdist_dendrogram`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_dendrogram.md),
[`plot_distance_distribution`](https://shihanli92.github.io/tcrdistR/reference/plot_distance_distribution.md),
[`plot_cdr3_logo`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md),
[`plot_junction_bars`](https://shihanli92.github.io/tcrdistR/reference/plot_junction_bars.md),
[`plot_gene_usage`](https://shihanli92.github.io/tcrdistR/reference/plot_gene_usage.md),
[`plot_tcr_scatter`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_scatter.md)

**I/O:**
[`read_tcr_table`](https://shihanli92.github.io/tcrdistR/reference/read_tcr_table.md),
[`read_airr`](https://shihanli92.github.io/tcrdistR/reference/read_airr.md),
[`read_adaptive`](https://shihanli92.github.io/tcrdistR/reference/read_adaptive.md),
[`read_10x`](https://shihanli92.github.io/tcrdistR/reference/read_10x.md)

## See also

[`vignette("tcrdistR-getting-started")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-getting-started.md)
for a tutorial introduction.

## Author

**Maintainer**: Shihan Li <shihanli1992@gmail.com>
