# tcrdistR 0.1.0

Initial release of tcrdistR, a C++-accelerated R package for TCR distance
computation and repertoire analysis targeting feature parity with Python
tcrdist3.

## TCRrep Object

* S4 `TCRrep` container with automatic deduplication and distance slots.

## Distance Computation

* Dense pairwise matrix (`tcrdist_matrix()`), sparse (`tcrdist_sparse()`),
  and rectangular (`tcrdist_rect()`) distance computation.
* Per-component decomposition via `components` parameter (CDR3, V-region,
  single-chain, custom combinations).
* Single-chain auto-detection for beta-only or alpha-only data.
* All core computation in C++ via Rcpp.

## Neighbor Search

* K-nearest neighbors (`tcrdist_knn()`), radius neighbors
  (`tcrdist_radius_neighbors()`).
* KNN from precomputed matrices (`knn_from_matrix()`) and PCA embeddings
  (`knn_from_pca()`).

## Dimensionality Reduction

* Kernel PCA with linear and Gaussian kernels
  (`compute_tcrdist_kernel_pca()`).
* UMAP via KNN or PCA path (`compute_tcrdist_umap()`).

## Statistical Analysis

* Neighborhood enrichment testing (`neighborhood_test()`).
* TCR clumping detection with background resampling (`find_clumping()`).
* Meta-clonotype discovery across subjects (`find_meta_clonotypes()`).
* Database matching with significance testing (`match_tcrs_to_db()`).

## Diversity Metrics

* Generalized Simpson's diversity (`tcr_diversity()`), fuzzy diversity
  (`tcr_fuzzy_diversity()`), richness (`tcr_richness()`), clonality
  (`tcr_clonality()`), Shannon entropy (`tcr_shannon_entropy()`), Gini
  index (`tcr_gini()`), repertoire overlap (`tcr_repertoire_overlap()`).

## Clustering

* Hierarchical clustering (`tcrdist_hclust()`) with k-cut and h-cut
  extraction (`cluster_tcrs()`).

## Visualization

* Heatmap, dendrogram, distance distribution, scatter plots, gene usage,
  CDR3 logos, junction bars, V/J gene logos, logo panels, network plots,
  CDR3 length distributions.

## I/O

* Readers for 10X, AIRR, Adaptive formats and auto-detection
  (`read_tcr_table()`).
* Format conversion from scRepertoire, scirpy, dandelion, tcrdist3
  (`as_tcr_df()`).

## Built-in Datasets

* `dash`: 1924 paired mouse TCRs (Dash et al., 2017).
* `flu`: Influenza-specific human TCR dataset.
