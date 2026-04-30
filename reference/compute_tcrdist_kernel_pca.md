# Kernel PCA on TCRdist distances

Computes a low-dimensional embedding of TCRs by applying kernel PCA to
the pairwise TCRdist distance matrix.

## Usage

``` r
compute_tcrdist_kernel_pca(
  tcr_df = NULL,
  organism = NULL,
  n_components = 50L,
  kernel = NULL,
  gaussian_kernel_sdev = 100,
  force_Dmax = NULL,
  method = c("auto", "eigen", "RSpectra"),
  dist_matrix = NULL
)
```

## Arguments

- tcr_df:

  A `data.frame` with at least columns `va`, `cdr3a`, `vb`, `cdr3b`.
  Optional if `dist_matrix` is provided.

- organism:

  Character string. Organism key, e.g. `"human"` or `"mouse"`. Optional
  if `dist_matrix` is provided.

- n_components:

  Integer. Maximum number of PCA components to return. Clamped to
  `nrow(tcr_df)`. Default `50L`.

- kernel:

  `NULL` (default linear kernel) or `"gaussian"`.

- gaussian_kernel_sdev:

  Numeric. Standard deviation parameter for the Gaussian kernel. Ignored
  unless `kernel = "gaussian"`. Default `100`.

- force_Dmax:

  Numeric or `NULL`. If non-`NULL`, use this value instead of `max(D)`
  when computing the default kernel. Ignored when `kernel = "gaussian"`.

- method:

  Character. Eigen-decomposition method: `"auto"` (default, uses
  [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html)
  when available for partial decomposition, falling back to
  [`base::eigen()`](https://rdrr.io/r/base/eigen.html)), `"eigen"`
  (always uses [`base::eigen()`](https://rdrr.io/r/base/eigen.html),
  same LAPACK as scipy.linalg.eigh), or `"RSpectra"` (always uses
  [`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html),
  same ARPACK as scipy.sparse.linalg.eigsh).

- dist_matrix:

  Optional precomputed distance matrix. If provided, `tcr_df` and
  `organism` are not used for distance computation.

## Value

A named list with elements:

- `embeddings`:

  Numeric matrix of dimensions N x n_components.

- `eigenvalues`:

  Numeric vector of retained positive eigenvalues (decreasing order).

- `n_components`:

  Integer. Number of components actually returned.

## Details

This implementation matches `scipy.linalg.eigh` (via sklearn's
`KernelPCA(kernel='precomputed')`). Both R's
[`base::eigen()`](https://rdrr.io/r/base/eigen.html) and scipy use the
same LAPACK `dsyevr` routine, so results are numerically identical to
~1e-10 tolerance.

Two kernel choices are supported:

- Default (`kernel = NULL`):

  Linear kernel: `gram = pmax(0, 1 - D / Dmax)` where
  `Dmax = force_Dmax %||% max(D)`.

- Gaussian (`kernel = "gaussian"`):

  RBF kernel: `gram = exp(-0.5 * (D / sdev)^2)`.

## See also

[`plot_tcr_scatter`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_scatter.md),
[`knn_from_pca`](https://shihanli92.github.io/tcrdistR/reference/knn_from_pca.md),
[`tcrdist_matrix`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_matrix.md)

## Examples

``` r
# \donttest{
tcrs <- data.frame(
    va    = c("TRAV1-1*01", "TRAV1-2*01", "TRAV1-1*01"),
    cdr3a = c("CAVRDSSYKLIF", "CAVRDSNYQLIW", "CAVRDSSYKLIF"),
    vb    = c("TRBV19*01", "TRBV28*01", "TRBV19*01"),
    cdr3b = c("CASSIRSSYEQYF", "CASSLGQAYEQYF", "CASSIRSYEQYF"),
    stringsAsFactors = FALSE
)
result <- compute_tcrdist_kernel_pca(tcrs, "human", n_components = 2L)
str(result)
#> List of 3
#>  $ embeddings  : num [1:3, 1:2] 0.47255 -0.92218 0.44963 0.23793 0.00398 ...
#>  $ eigenvalues : num [1:2] 1.276 0.115
#>  $ n_components: int 2
# }
```
