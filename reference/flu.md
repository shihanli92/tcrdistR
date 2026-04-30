# Human influenza TCR dataset from VDJdb

A dataset of 2271 paired alpha-beta human T-cell receptors specific to 4
influenza epitopes, compiled from the VDJdb database (September 2020
release). This dataset covers both MHC class I and class II restricted
responses from 63 subjects and is useful for benchmarking human TCR
distance calculations.

## Usage

``` r
flu
```

## Format

A data.frame with 2271 rows and 13 columns:

- subject:

  Subject identifier from the original study (may be `NA` if not
  annotated).

- epitope:

  Epitope peptide sequence (e.g., `"GILGFVFTL"`).

- epitope_gene:

  Source gene of the epitope (e.g., `"M1"`).

- mhc_a:

  MHC alpha chain allele (e.g., `"HLA-A*02"`).

- mhc_b:

  MHC beta chain allele (e.g., `"B2M"`).

- mhc_class:

  MHC class: `"MHCI"` or `"MHCII"`.

- count:

  Clone count (set to 1 for all entries).

- va:

  V-alpha gene with allele (e.g., `"TRAV12-2*01"`).

- ja:

  J-alpha gene with allele (e.g., `"TRAJ33*01"`).

- cdr3a:

  CDR3-alpha amino acid sequence.

- vb:

  V-beta gene with allele (e.g., `"TRBV19*01"`).

- jb:

  J-beta gene with allele (e.g., `"TRBJ2-7*01"`).

- cdr3b:

  CDR3-beta amino acid sequence.

## Source

VDJdb: <https://vdjdb.cdr3.net/>

Bagaev et al. (2020). VDJdb in 2019: database extension, new analysis
infrastructure and a T-cell receptor motif compendium. *Nucleic Acids
Research*, 48(D1), D1057–D1062.
[doi:10.1093/nar/gkz874](https://doi.org/10.1093/nar/gkz874)

## Examples

``` r
data(flu)
dim(flu)         # 2271 x 13
#> [1] 2271   13
table(flu$epitope)
#> 
#>     CVNGSCFTV  DATYQRTRALVR     GILGFVFTL PKYVKQNTLKLAT 
#>            14           109          2078            70 
```
