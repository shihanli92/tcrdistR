# Single-linkage clustering of clumped TCRs

Given a set of clumped TCR indices and their neighbor sets, performs
single-linkage clustering by iteratively propagating the smallest
neighbor index until convergence.

## Usage

``` r
.single_linkage_clumping(all_clumped_nbrs, num_clones, is_clumped)
```

## Arguments

- all_clumped_nbrs:

  Named list (character keys are 0-based clone indices). Each value is
  an integer vector of 0-based neighbor indices.

- num_clones:

  Integer. Total number of clones.

- is_clumped:

  Logical vector of length `num_clones`.

## Value

Integer vector of length `num_clones`. 0 means not clumped; positive
integers are clumping group IDs (reordered by decreasing group size).
