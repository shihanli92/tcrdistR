# CD8 logistic regression score for paired TCRs

Scores each TCR using a logistic regression model trained on sorted
CD4/CD8 bulk TCR data. The alpha-chain and beta-chain models are fit
independently and their scores are averaged (each weighted 0.5).

## Usage

``` r
make_cd8_score_table_column(tcr_df, use_sigmoid = FALSE)
```

## Arguments

- tcr_df:

  A data.frame with columns `va`, `ja`, `cdr3a`, `vb`, `jb`, `cdr3b`.
  Gene names should include allele suffixes (e.g. `"TRAV1-2*01"`).

- use_sigmoid:

  Logical. If `TRUE`, apply the sigmoid function \\1 / (1 + exp(-x))\\
  to each chain's score before averaging. Default is `FALSE`.

## Value

A numeric vector of CD8 logistic regression scores, one per row of
`tcr_df`.

## Details

The model encodes V-gene, J-gene, CDR3 length, and positional amino acid
features using one-hot and k-hot encodings, then computes a dot product
with pre-trained weights.

## See also

[`match_tcrs_to_db`](https://shihanli92.github.io/tcrdistR/reference/match_tcrs_to_db.md)

## Examples

``` r
# \donttest{
tcr_df <- data.frame(
    va = "TRAV1-2*01", ja = "TRAJ33*01",
    cdr3a = "CAVMDSSYKLIF",
    vb = "TRBV6-4*01", jb = "TRBJ2-1*01",
    cdr3b = "CASSLAPGATNEKLFF",
    stringsAsFactors = FALSE
)
scores <- make_cd8_score_table_column(tcr_df)
# }
```
