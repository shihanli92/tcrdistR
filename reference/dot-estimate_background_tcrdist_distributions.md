# Estimate background paired TCRdist distributions

For each foreground TCR, estimates the probability of seeing a paired
TCRdist score at or below each integer distance from 0 to `max_dist`
under a null model where alpha and beta chains are independently drawn
from shuffled backgrounds.

## Usage

``` r
.estimate_background_tcrdist_distributions(
  organism,
  tcr_df,
  max_dist,
  num_random_samples = 50000L,
  pseudocount = 0.25,
  preserve_vj_pairings = FALSE,
  bg_tcrs = NULL
)
```

## Arguments

- organism:

  Character string. Organism key (e.g. `"human"`).

- tcr_df:

  A data.frame with columns
  `va, ja, cdr3a, cdr3a_nucseq, vb, jb, cdr3b, cdr3b_nucseq`.

- max_dist:

  Integer. Maximum paired TCRdist to consider.

- num_random_samples:

  Integer. Number of random background chains to generate per chain
  type.

- pseudocount:

  Numeric. Pseudocount added before normalizing.

- preserve_vj_pairings:

  Logical. Preserve V-J pairings in resampling.

- bg_tcrs:

  Optional data.frame used for background generation.

## Value

Numeric matrix of dimensions `nrow(tcr_df) x (max_dist + 1)`.
