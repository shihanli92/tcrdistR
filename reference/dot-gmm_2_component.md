# Fit a 2-component Gaussian mixture via EM and find the crossover point

Runs expectation-maximization for a mixture of two univariate Gaussians.
Returns the point where the posterior probabilities are equal (the
"crossover"), which corresponds to the valley between components.

## Usage

``` r
.gmm_2_component(x, max_iter = 200L, tol = 1e-06)
```

## Arguments

- x:

  Numeric vector of distances.

- max_iter:

  Integer. Maximum EM iterations.

- tol:

  Numeric. Convergence tolerance on log-likelihood.

## Value

Numeric scalar (crossover point) or `NULL` on failure.
