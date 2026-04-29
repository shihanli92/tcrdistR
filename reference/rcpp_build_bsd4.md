# Build the BSD4 amino acid distance matrix (C++ constexpr values)

Returns the 20x20 BSD4 distance matrix derived from BLOSUM62 using the
CoNGA transformation:

- diagonal (a == b): 0

- off-diagonal, BLOSUM62 \< 0: 4

- off-diagonal, BLOSUM62 \>= 0: 4 - BLOSUM62

Row and column names follow alphabetical single-letter amino acid order
(A C D E F G H I K L M N P Q R S T V W Y), matching `AMINO_ACIDS` in
rconga.

## Usage

``` r
rcpp_build_bsd4()
```

## Value

A named 20x20 `NumericMatrix` of BSD4 distances.

## Details

The values are compiled into the package as `constexpr` C++ arrays; this
function exposes them to R for verification and for use by any R code
that needs the matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
  bsd4 <- rcpp_build_bsd4()
  bsd4["A", "S"]  # 3 (BLOSUM62[A,S]=1, BSD4=4-1=3)
  bsd4["A", "A"]  # 0 (diagonal)
  bsd4["A", "C"]  # 4 (BLOSUM62=-something negative)
} # }
```
