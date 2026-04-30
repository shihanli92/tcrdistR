# tcrdistR

C++-accelerated TCR distance calculations for T-cell receptor repertoire
analysis. Computes pairwise TCRdist distances incorporating V-region and
CDR3 sequence comparisons using BLOSUM62-derived substitution matrices.
Supports paired alpha-beta and single-chain (beta-only or alpha-only)
input, per-component distance decomposition, and targets feature parity
with Python [tcrdist3](https://github.com/kmayerb/tcrdist3).

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("shihanli92/tcrdistR")
```

## Quick Start

``` r
library(tcrdistR)
data(dash)

# Build a TCRrep object (deduplicates identical clones per subject)
rep <- TCRrep(dash, organism = "mouse", compute_distances = TRUE)
rep
#> TCRrep with 1888 clones (mouse, AB chains)
#>   paired_dist: 1888 x 1888

# Pairwise distance heatmap
plot_tcrdist_heatmap(rep@paired_dist[1:40, 1:40])

# Per-component distances (CDR3-only, V-region-only, etc.)
d_cdr3 <- tcrdist_matrix(rep@clone_df, "mouse", components = "cdr3")
d_v    <- tcrdist_matrix(rep@clone_df, "mouse", components = "v_region")

# Single-chain mode (beta-only input, no alpha columns needed)
beta_only <- rep@clone_df[, c("vb", "cdr3b")]
d_beta <- tcrdist_matrix(beta_only, "mouse")

# Kernel PCA
pca <- compute_tcrdist_kernel_pca(
  rep@clone_df, rep@organism, n_components = 50L
)
plot_tcr_scatter(
  pca$embeddings[, 1:2],
  color_by = rep@clone_df$epitope,
  title = "Kernel PCA",
  point_size = 1.5
)

# UMAP (uses TCRdist KNN neighbors directly)
umap <- compute_tcrdist_umap(rep@clone_df, rep@organism, seed = 42)
plot_tcr_scatter(
  umap$embeddings,
  color_by = rep@clone_df$epitope,
  axis_label_prefix = "UMAP",
  point_size = 1.5
)
```

## Documentation

Full reference and vignettes at
**<https://shihanli92.github.io/tcrdistR/>**.

- [`vignette("tcrdistR-getting-started")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-getting-started.md)
  – Installation and core functions
- [`vignette("tcrdistR-tcrrep-workflow")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-tcrrep-workflow.md)
  – End-to-end analysis with the TCRrep object
- [`vignette("tcrdistR-advanced")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-advanced.md)
  – Clumping, meta-clonotypes, clustering, DB matching
- [`vignette("tcrdistR-visualization")`](https://shihanli92.github.io/tcrdistR/articles/tcrdistR-visualization.md)
  – Visualization guide

## License

MIT
