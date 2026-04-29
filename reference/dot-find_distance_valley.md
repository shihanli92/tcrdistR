# Find the valley between the two major peaks of a bimodal distribution

Uses kernel density estimation to locate the minimum between the two
dominant peaks of the typically bimodal TCRdist distribution. First
identifies the two tallest local maxima, then finds the deepest local
minimum between them. Falls back to a 2-component Gaussian mixture model
(EM) when KDE peak detection fails.

## Usage

``` r
.find_distance_valley(distances, n_points = 512L, adjust = 1)
```

## Arguments

- distances:

  Numeric vector of pairwise distances.

- n_points:

  Integer. Number of density estimation points.

- adjust:

  Numeric. Bandwidth adjustment for
  [`density()`](https://rdrr.io/r/stats/density.html).

## Value

Numeric scalar. The distance at the valley.
