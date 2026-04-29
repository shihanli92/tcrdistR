# Run graph-based clustering on a fuzzy simplicial set graph

Performs Leiden or Louvain community detection on the KNN graph. When
`method` is `NULL`, Leiden is tried first with Louvain as fallback. Port
of rconga's `.run_clustering()`.

## Usage

``` r
.run_clustering(knn_graph, resolution = 1, method = NULL)
```

## Arguments

- knn_graph:

  Sparse matrix (N x N) from
  [`.build_knn_graph()`](https://shihanli92.github.io/tcrdistR/reference/dot-build_knn_graph.md).

- resolution:

  Numeric. Resolution parameter. Default `1.0`.

- method:

  Character or `NULL`. `"leiden"`, `"louvain"`, or `NULL` (try Leiden
  first).

## Value

Integer vector of length N with 0-based cluster IDs.
