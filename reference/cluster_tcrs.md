# Cluster TCRs using various algorithms

Assigns TCR clonotypes to clusters using one of five methods:
hierarchical clustering, Leiden or Louvain community detection, DBSCAN
density-based clustering, or k-medoids (PAM).

## Usage

``` r
cluster_tcrs(
  tcr_df = NULL,
  organism = NULL,
  method = c("hierarchical", "leiden", "louvain", "dbscan", "kmedoids"),
  dist_matrix = NULL,
  k = NULL,
  h = NULL,
  hclust_method = "average",
  resolution = 1,
  n_neighbors = 10L,
  eps = NULL,
  min_pts = 5L
)
```

## Arguments

- tcr_df:

  Data.frame with TCR columns. Required unless `dist_matrix` is
  provided.

- organism:

  Character string (`"human"` or `"mouse"`). Required unless
  `dist_matrix` is provided.

- method:

  Clustering method. One of `"hierarchical"` (default), `"leiden"`,
  `"louvain"`, `"dbscan"`, or `"kmedoids"`.

- dist_matrix:

  Optional precomputed N x N distance matrix. If provided, `tcr_df` and
  `organism` are only needed for graph-based methods (Leiden/Louvain)
  when KNN must be computed.

- k:

  Integer. Number of clusters for `"hierarchical"` and `"kmedoids"`. For
  hierarchical, exactly one of `k` or `h` must be specified. For
  kmedoids, required.

- h:

  Numeric. Height for dendrogram cutting (hierarchical only).

- hclust_method:

  Character. Agglomeration method for
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html). Default
  `"average"` (UPGMA).

- resolution:

  Numeric. Resolution parameter for Leiden/Louvain. Higher values yield
  more clusters. Default `1.0`.

- n_neighbors:

  Integer. Number of nearest neighbors for KNN graph construction
  (Leiden/Louvain). Default `10L`.

- eps:

  Numeric or `NULL`. DBSCAN neighborhood radius. If `NULL`,
  auto-detected from the k-distance knee.

- min_pts:

  Integer. DBSCAN minimum points for core point. Default `5L`.

## Value

An integer vector of cluster assignments (length N).

- For `hierarchical`, `leiden`, `louvain`, and `kmedoids`: 1-based
  cluster IDs.

- For `dbscan`: 0-based (0 = noise/unassigned), with cluster IDs
  starting at 1.

## Details

Distances can be precomputed via `dist_matrix` or computed internally
from `tcr_df` and `organism`.

## See also

[`tcrdist_hclust`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_hclust.md),
[`neighborhood_test`](https://shihanli92.github.io/tcrdistR/reference/neighborhood_test.md),
[`compute_tcrdist_umap`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_umap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(dash)
sub <- dash[1:50, ]

# Hierarchical (default)
clusters <- cluster_tcrs(sub, "mouse", k = 5)

# K-medoids
km <- cluster_tcrs(sub, "mouse", method = "kmedoids", k = 4)

# DBSCAN with auto-detected eps
db <- cluster_tcrs(sub, "mouse", method = "dbscan")

# Leiden on KNN graph
lei <- cluster_tcrs(sub, "mouse", method = "leiden", resolution = 1.0)

# Precomputed distance matrix
dm <- tcrdist_matrix(sub, "mouse")
clusters <- cluster_tcrs(dist_matrix = dm, method = "kmedoids", k = 3)
} # }
```
