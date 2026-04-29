# Compute UMAP embedding from TCRdist kernel PCA

Reduces kernel PCA embeddings (or raw TCR data) to a low-dimensional
UMAP representation suitable for visualization with
[`plot_tcr_scatter`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_scatter.md).
Uses the uwot package for the UMAP computation.

## Usage

``` r
compute_tcrdist_umap(
  tcr_df = NULL,
  organism = NULL,
  pca_embeddings = NULL,
  n_components_pca = 50L,
  n_components = 2L,
  n_neighbors = 15L,
  min_dist = 0.1,
  metric = "euclidean",
  n_threads = 1L,
  seed = NULL,
  ...
)
```

## Arguments

- tcr_df:

  Data.frame with TCR columns (`va`, `cdr3a`, `vb`, `cdr3b`). Used only
  if `pca_embeddings` is `NULL`.

- organism:

  Character string (`"human"` or `"mouse"`). Required when `tcr_df` is
  used.

- pca_embeddings:

  Numeric matrix (N x D). Pre-computed kernel PCA embeddings. If
  provided, `tcr_df` and `organism` are ignored.

- n_components_pca:

  Integer. Number of kernel PCA components to compute when using the
  `tcr_df` input path. Default `50L`.

- n_components:

  Integer. Number of UMAP output dimensions. Default `2L`.

- n_neighbors:

  Integer. Size of local neighborhood for UMAP manifold approximation.
  Default `15L`.

- min_dist:

  Numeric. Minimum distance between embedded points. Controls how
  tightly UMAP packs points together. Default `0.1`.

- metric:

  Character string. Distance metric for UMAP neighbor search. Default
  `"euclidean"`.

- n_threads:

  Integer. Number of threads for neighbor search and optimization.
  Default `1L`.

- seed:

  Integer or `NULL`. Random seed for reproducibility. If non-`NULL`,
  [`set.seed()`](https://rdrr.io/r/base/Random.html) is called before
  UMAP computation. Default `NULL`.

- ...:

  Additional arguments passed to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

## Value

A named list:

- `embeddings`:

  Numeric matrix (N x `n_components`). UMAP coordinates.

- `pca_embeddings`:

  Numeric matrix. The PCA input used.

- `n_components`:

  Integer. Number of UMAP dimensions.

## Details

Two input modes are supported:

1.  **From pre-computed PCA**: supply `pca_embeddings` (an N x D matrix,
    e.g. from `compute_tcrdist_kernel_pca()$embeddings`).

2.  **From raw TCR data**: supply `tcr_df` and `organism`. Kernel PCA is
    computed internally with `n_components_pca` components before UMAP.

## See also

[`compute_tcrdist_kernel_pca`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_kernel_pca.md),
[`plot_tcr_scatter`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_scatter.md)

## Examples

``` r
# \donttest{
data(dash)
# From raw data (computes kernel PCA internally)
umap <- compute_tcrdist_umap(dash[1:100, ], "mouse", seed = 42)
dim(umap$embeddings)  # 100 x 2
#> [1] 100   2

# From pre-computed PCA
pca <- compute_tcrdist_kernel_pca(dash[1:100, ], "mouse", n_components = 20L)
umap <- compute_tcrdist_umap(pca_embeddings = pca$embeddings, seed = 42)
# }
```
