## Create the `dash` package dataset from the DASH fixture CSV.
##
## Source: Dash et al. (2017), "Quantifiable predictive features define
## epitope-specific T cell receptor repertoires", Nature 547, 89-93.
## doi:10.1038/nature22383
##
## The original CSV lives at tests/testthat/fixtures/dash.csv and contains
## 1924 paired alpha-beta mouse TCRs across 7 epitopes from 78 subjects.

raw <- read.csv(
  "tests/testthat/fixtures/dash.csv",
  stringsAsFactors = FALSE
)

dash <- data.frame(
  subject        = raw$subject,
  epitope        = raw$epitope,
  count          = raw$count,
  va             = raw$v_a_gene,
  ja             = raw$j_a_gene,
  cdr3a          = raw$cdr3_a_aa,
  cdr3a_nucseq   = raw$cdr3_a_nucseq,
  vb             = raw$v_b_gene,
  jb             = raw$j_b_gene,
  cdr3b          = raw$cdr3_b_aa,
  cdr3b_nucseq   = raw$cdr3_b_nucseq,
  clone_id       = raw$clone_id,
  stringsAsFactors = FALSE
)

usethis::use_data(dash, overwrite = TRUE, compress = "xz")
