# Standardize a TCR data.frame to tcrdistR column names

Converts a data.frame from common TCR analysis tools to tcrdistR's
canonical column names (`va`, `cdr3a`, `vb`, `cdr3b`). Supports
auto-detection of column naming patterns from scRepertoire, scirpy,
dandelion, and tcrdist3, as well as custom column mappings.

## Usage

``` r
as_tcr_df(
  df,
  col_map = NULL,
  format = NULL,
  normalize_genes = TRUE,
  drop_incomplete = TRUE
)
```

## Arguments

- df:

  A `data.frame` containing TCR data.

- col_map:

  Named character vector mapping tcrdistR canonical names to the column
  names in `df`. For example,
  `c(va = "alpha_v_gene", cdr3a = "alpha_cdr3")`. If `NULL`, the format
  is auto-detected.

- format:

  Character string or `NULL`. Force a specific format instead of
  auto-detecting:

  `"screpertoire"`

  :   scRepertoire `CTgene`/`CTaa` concatenated columns.

  `"scirpy"`, `"dandelion"`

  :   scirpy/dandelion `IR_VJ_1_*`/`IR_VDJ_1_*` columns.

  `"tcrdist3"`

  :   tcrdist3 `v_a_gene`/`cdr3_a_aa` columns.

  If `NULL` (default), auto-detection is used.

- normalize_genes:

  Logical. If `TRUE` (default), gene names without an allele suffix get
  `*01` appended.

- drop_incomplete:

  Logical. If `TRUE` (default), rows with `NA` or empty strings in
  required chain columns (`va`, `cdr3a`, `vb`, `cdr3b`) are dropped with
  a message.

## Value

A `data.frame` with tcrdistR canonical column names. Extra columns from
the input are preserved.

## See also

[`read_tcr_table`](https://shihanli92.github.io/tcrdistR/reference/read_tcr_table.md),
[`TCRrep`](https://shihanli92.github.io/tcrdistR/reference/TCRrep.md)

## Examples

``` r
# Custom column mapping
df <- data.frame(
  alpha_v = "TRAV1-1", alpha_cdr3 = "CAVRDSSYKLIF",
  beta_v  = "TRBV5-1", beta_cdr3  = "CASSIRSSYEQYF"
)
as_tcr_df(df, col_map = c(va = "alpha_v", cdr3a = "alpha_cdr3",
                           vb = "beta_v",  cdr3b = "beta_cdr3"))
#> as_tcr_df: detected format 'custom'
#>           va        cdr3a         vb         cdr3b
#> 1 TRAV1-1*01 CAVRDSSYKLIF TRBV5-1*01 CASSIRSSYEQYF

# scirpy / dandelion format (auto-detected)
df <- data.frame(
  IR_VJ_1_v_call = "TRAV1-1", IR_VJ_1_junction_aa = "CAVRDSSYKLIF",
  IR_VDJ_1_v_call = "TRBV5-1", IR_VDJ_1_junction_aa = "CASSIRSSYEQYF"
)
as_tcr_df(df)
#> as_tcr_df: detected format 'scirpy'
#>           va        cdr3a         vb         cdr3b
#> 1 TRAV1-1*01 CAVRDSSYKLIF TRBV5-1*01 CASSIRSSYEQYF

# scRepertoire format (auto-detected)
df <- data.frame(
  CTgene = "TRAV1-1.TRAJ33.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
  CTaa   = "CAVRDSSYKLIF_CASSIRSSYEQYF"
)
as_tcr_df(df)
#> as_tcr_df: detected format 'screpertoire'
#>           va        ja        cdr3a         vb         jb         cdr3b
#> 1 TRAV1-1*01 TRAJ33*01 CAVRDSSYKLIF TRBV5-1*01 TRBJ2-7*01 CASSIRSSYEQYF
```
