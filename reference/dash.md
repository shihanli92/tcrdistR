# DASH dataset: paired alpha-beta mouse TCRs across 7 epitopes

A dataset of 1924 paired alpha-beta T-cell receptors from mice
responding to 7 viral epitopes, collected from 78 subjects. This is the
benchmark dataset from Dash et al. (2017) and is widely used for
evaluating TCR distance metrics.

## Usage

``` r
dash
```

## Format

A data.frame with 1924 rows and 12 columns:

- subject:

  Subject identifier (e.g., `"mouse_subject0050"`).

- epitope:

  Epitope specificity: F2, M38, M45, m139, NP, PA, or PB1.

- count:

  Clone count (number of cells observed for this clonotype).

- va:

  V-alpha gene with allele (e.g., `"TRAV7-3*01"`).

- ja:

  J-alpha gene with allele (e.g., `"TRAJ33*01"`).

- cdr3a:

  CDR3-alpha amino acid sequence (e.g., `"CAVSLDSNYQLIW"`).

- cdr3a_nucseq:

  CDR3-alpha nucleotide sequence.

- vb:

  V-beta gene with allele (e.g., `"TRBV13-1*01"`).

- jb:

  J-beta gene with allele (e.g., `"TRBJ2-3*01"`).

- cdr3b:

  CDR3-beta amino acid sequence (e.g., `"CASSDFDWGGDAETLYF"`).

- cdr3b_nucseq:

  CDR3-beta nucleotide sequence.

- clone_id:

  Unique clone identifier.

## Source

Dash et al. (2017). Quantifiable predictive features define
epitope-specific T cell receptor repertoires. *Nature*, 547, 89–93.
[doi:10.1038/nature22383](https://doi.org/10.1038/nature22383)

## Examples

``` r
data(dash)
dim(dash)       # 1924 x 12
#> [1] 1924   12
table(dash$epitope)
#> 
#>   F2  M38  M45   NP   PA  PB1 m139 
#>  117  158  291  305  324  642   87 
```
