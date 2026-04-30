# Changelog

## tcrdistR 0.1.0

Initial release of tcrdistR, a C++-accelerated R package for TCR
distance computation and repertoire analysis targeting feature parity
with Python tcrdist3.

### TCRrep Object

- S4 `TCRrep` container with automatic deduplication and distance slots.

### Distance Computation

- Dense pairwise matrix
  ([`tcrdist_matrix()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)),
  sparse
  ([`tcrdist_sparse()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_sparse.md)),
  and rectangular
  ([`tcrdist_rect()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_rect.md))
  distance computation.
- Per-component decomposition via `components` parameter (CDR3,
  V-region, single-chain, custom combinations).
- Single-chain auto-detection for beta-only or alpha-only data.
- All core computation in C++ via Rcpp.

### Neighbor Search

- K-nearest neighbors
  ([`tcrdist_knn()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_knn.md)),
  radius neighbors
  ([`tcrdist_radius_neighbors()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_radius_neighbors.md)).
- KNN from precomputed matrices
  ([`knn_from_matrix()`](https://shihanli92.github.io/tcrdistR/reference/knn_from_matrix.md))
  and PCA embeddings
  ([`knn_from_pca()`](https://shihanli92.github.io/tcrdistR/reference/knn_from_pca.md)).

### Dimensionality Reduction

- Kernel PCA with linear and Gaussian kernels
  ([`compute_tcrdist_kernel_pca()`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_kernel_pca.md)).
- UMAP via KNN or PCA path
  ([`compute_tcrdist_umap()`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_umap.md)).

### Statistical Analysis

- Neighborhood enrichment testing
  ([`neighborhood_test()`](https://shihanli92.github.io/tcrdistR/reference/neighborhood_test.md)).
- TCR clumping detection with background resampling
  ([`find_clumping()`](https://shihanli92.github.io/tcrdistR/reference/find_clumping.md)).
- Meta-clonotype discovery across subjects
  ([`find_meta_clonotypes()`](https://shihanli92.github.io/tcrdistR/reference/find_meta_clonotypes.md)).
- Database matching with significance testing
  ([`match_tcrs_to_db()`](https://shihanli92.github.io/tcrdistR/reference/match_tcrs_to_db.md)).

### Diversity Metrics

- Generalized Simpson’s diversity
  ([`tcr_diversity()`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md)),
  fuzzy diversity
  ([`tcr_fuzzy_diversity()`](https://shihanli92.github.io/tcrdistR/reference/tcr_fuzzy_diversity.md)),
  richness
  ([`tcr_richness()`](https://shihanli92.github.io/tcrdistR/reference/tcr_richness.md)),
  clonality
  ([`tcr_clonality()`](https://shihanli92.github.io/tcrdistR/reference/tcr_clonality.md)),
  Shannon entropy
  ([`tcr_shannon_entropy()`](https://shihanli92.github.io/tcrdistR/reference/tcr_shannon_entropy.md)),
  Gini index
  ([`tcr_gini()`](https://shihanli92.github.io/tcrdistR/reference/tcr_gini.md)),
  repertoire overlap
  ([`tcr_repertoire_overlap()`](https://shihanli92.github.io/tcrdistR/reference/tcr_repertoire_overlap.md)).

### Clustering

- Hierarchical clustering
  ([`tcrdist_hclust()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_hclust.md))
  with k-cut and h-cut extraction
  ([`cluster_tcrs()`](https://shihanli92.github.io/tcrdistR/reference/cluster_tcrs.md)).

### Visualization

- Heatmap, dendrogram, distance distribution, scatter plots, gene usage,
  CDR3 logos, junction bars, V/J gene logos, logo panels, network plots,
  CDR3 length distributions.

### I/O

- Readers for 10X, AIRR, Adaptive formats and auto-detection
  ([`read_tcr_table()`](https://shihanli92.github.io/tcrdistR/reference/read_tcr_table.md)).
- Format conversion from scRepertoire, scirpy, dandelion, tcrdist3
  ([`as_tcr_df()`](https://shihanli92.github.io/tcrdistR/reference/as_tcr_df.md)).

### Built-in Datasets

- `dash`: 1924 paired mouse TCRs (Dash et al., 2017).
- `flu`: Influenza-specific human TCR dataset.
