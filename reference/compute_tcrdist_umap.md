# Compute UMAP embedding from TCRdist data

Two input modes are supported:

## Usage

``` r
compute_tcrdist_umap(
  tcr_df = NULL,
  organism = NULL,
  pca_embeddings = NULL,
  n_components = 2L,
  n_neighbors = 15L,
  min_dist = 0.1,
  spread = 1,
  metric = "euclidean",
  n_threads = 1L,
  seed = NULL,
  cluster = FALSE,
  clustering_resolution = 1,
  clustering_method = NULL,
  ...
)
```

## Arguments

- tcr_df:

  Data.frame with TCR columns (`va`, `cdr3a`, `vb`, `cdr3b`; plus `ja`,
  `jb` for group masking).

- organism:

  Character string (`"human"` or `"mouse"`). Required when `tcr_df` is
  used.

- pca_embeddings:

  Numeric matrix (N x D). Pre-computed kernel PCA embeddings. If
  provided, the PCA path is used.

- n_components:

  Integer. Number of UMAP output dimensions. Default `2L`.

- n_neighbors:

  Integer. Number of nearest neighbors. Default `15L`.

- min_dist:

  Numeric. UMAP min_dist. Default `0.1`.

- spread:

  Numeric. UMAP spread parameter. Default `1.0`.

- metric:

  Character. Distance metric for UMAP neighbor search (PCA path only).
  Default `"euclidean"`.

- n_threads:

  Integer. Threads for UMAP optimization. Default `1L`.

- seed:

  Integer or `NULL`. Random seed. Default `NULL`.

- cluster:

  Logical. If `TRUE`, run graph-based clustering on the KNN graph (KNN
  path only; requires igraph). Default `FALSE`.

- clustering_resolution:

  Numeric. Resolution for community detection. Default `1.0`.

- clustering_method:

  Character or `NULL`. `"leiden"`, `"louvain"`, or `NULL` (try Leiden
  first). Default `NULL`.

- ...:

  Additional arguments passed to
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.html) (PCA
  path only).

## Value

A named list. Elements depend on the input path:

**KNN path** (`tcr_df` + `organism`):

- `embeddings`:

  Numeric matrix (N x `n_components`).

- `knn_indices`:

  Integer matrix (N x K). 1-based.

- `knn_distances`:

  Numeric matrix (N x K).

- `knn_graph`:

  Sparse `dgCMatrix` (N x N). Fuzzy simplicial set.

- `clusters`:

  Integer vector (0-based) or `NULL`.

- `nndists`:

  Numeric vector. Weighted NN distances.

- `n_neighbors`:

  Final K (may have been increased for connectivity).

- `n_components`:

  Integer.

- `method`:

  `"knn"`.

**PCA path** (`pca_embeddings`):

- `embeddings`:

  Numeric matrix (N x `n_components`).

- `pca_embeddings`:

  The PCA input matrix.

- `n_components`:

  Integer.

- `method`:

  `"pca"`.

## Details

1.  **KNN path** (recommended): supply `tcr_df` and `organism`. Computes
    TCRdist K-nearest-neighbors with group masking (clones sharing an
    identical alpha or beta chain are excluded from each other's
    neighborhoods), builds a fuzzy simplicial set graph, and runs UMAP
    from the precomputed KNN via
    [`umap`](https://jlmelville.github.io/uwot/reference/umap.html).
    This preserves the TCRdist metric faithfully.

2.  **PCA path**: supply `pca_embeddings` (an N x D matrix, e.g. from
    `compute_tcrdist_kernel_pca()$embeddings`). Runs standard UMAP in
    Euclidean space on the PCA coordinates.

When `cluster = TRUE`, Leiden (or Louvain) community detection is
performed on the fuzzy KNN graph (KNN path only; requires igraph).

## See also

[`tcrdist_knn`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_knn.md),
[`setup_tcr_groups`](https://shihanli92.github.io/tcrdistR/reference/setup_tcr_groups.md),
[`compute_tcrdist_kernel_pca`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_kernel_pca.md),
[`plot_tcr_scatter`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_scatter.md)

## Examples

``` r
# \donttest{
data(dash)
sub <- dash[1:200, ]

# KNN path (recommended): uses TCRdist directly
umap <- compute_tcrdist_umap(sub, "mouse", seed = 42)
dim(umap$embeddings)  # 200 x 2
#> [1] 200   2

# With clustering
umap <- compute_tcrdist_umap(sub, "mouse", seed = 42, cluster = TRUE)
table(umap$clusters)
#> 
#>  0  1  2  3  4  5  6  7  8  9 10 
#> 27  9 22 35 11  9 21 28 14 16  8 

# PCA path: from pre-computed kernel PCA
pca <- compute_tcrdist_kernel_pca(sub, "mouse", n_components = 20L)
umap <- compute_tcrdist_umap(pca_embeddings = pca$embeddings, seed = 42)
# }
```
