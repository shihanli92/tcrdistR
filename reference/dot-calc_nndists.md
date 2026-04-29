# Weighted nearest-neighbor distances

Computes a weighted average of each row of the sorted KNN distance
matrix. Weights decrease linearly from 1.0 (nearest) to 1/K (farthest).
Port of rconga's `.calc_nndists()`.

## Usage

``` r
.calc_nndists(knn_distances)
```

## Arguments

- knn_distances:

  Numeric matrix (N x K). Sorted ascending.

## Value

Numeric vector of length N.
