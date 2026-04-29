# Double-center a kernel (Gram) matrix

Applies the standard double-centering transform used by kernel PCA
(matches `sklearn.preprocessing.KernelCenterer`): \$\$K_c = K - 1_n K -
K 1_n + 1_n K 1_n\$\$ where \\1_n\\ is an n-by-n matrix of \\1/n\\.

## Usage

``` r
.center_kernel_matrix(K)
```

## Arguments

- K:

  Numeric square matrix (the kernel / Gram matrix).

## Value

The double-centered kernel matrix (same dimensions as `K`).
