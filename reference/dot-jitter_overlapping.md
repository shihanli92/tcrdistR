# Jitter overlapping layout coordinates

Finds groups of nodes that share identical (x, y) positions and spreads
each group into a small circle. The radius is proportional to the
overall layout extent so the offset is visible but not disruptive.

## Usage

``` r
.jitter_overlapping(coords, frac = 0.015)
```

## Arguments

- coords:

  Numeric matrix (N x 2) of layout coordinates.

- frac:

  Numeric. Jitter radius as a fraction of the layout extent. Default
  `0.015`.

## Value

Numeric matrix (N x 2) with adjusted coordinates.
