## Create the `flu` package dataset from VDJdb influenza data.
##
## Source: VDJdb (https://vdjdb.cdr3.net/) — Influenza A/B TCR records,
## downloaded September 2020. Contains paired alpha-beta human TCRs with
## known epitope specificities.
##
## Original file: influenza_2020-SEP.tsv from
##   https://www.dropbox.com/s/mkjdaygdl41piw6/influenza_2020-SEP.tsv.zip?dl=1
##
## The VDJdb format stores alpha and beta chains as separate rows linked by
## complex.id. Rows with complex.id = 0 are unpaired and excluded.

raw <- read.delim(
  "data-raw/influenza_2020-SEP.tsv",
  stringsAsFactors = FALSE,
  quote = ""
)

# Keep only paired human entries (complex.id != 0, HomoSapiens only)
paired <- raw[raw$complex.id != 0 & raw$Species == "HomoSapiens", ]

# Split by chain
alpha <- paired[paired$Gene == "TRA", ]
beta  <- paired[paired$Gene == "TRB", ]

# Merge on complex.id
merged <- merge(alpha, beta, by = "complex.id", suffixes = c("_a", "_b"))

# Extract subject.id from the Meta JSON (simple regex — avoids jsonlite dep)
extract_meta_field <- function(meta, field) {
  pattern <- paste0('"', field, '":\\s*"([^"]*)"')
  m <- regmatches(meta, regexec(pattern, meta))
  vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
}

subject_id <- extract_meta_field(merged$Meta_b, "subject.id")
subject_id[subject_id == ""] <- NA_character_

# Build the output data.frame
flu <- data.frame(
  subject     = subject_id,
  epitope     = merged$Epitope_b,
  epitope_gene = merged$Epitope.gene_b,
  mhc_a       = merged$MHC.A_b,
  mhc_b       = merged$MHC.B_b,
  mhc_class   = merged$MHC.class_b,
  count       = 1L,
  va          = merged$V_a,
  ja          = merged$J_a,
  cdr3a       = merged$CDR3_a,
  vb          = merged$V_b,
  jb          = merged$J_b,
  cdr3b       = merged$CDR3_b,
  stringsAsFactors = FALSE
)

# Sort by epitope then subject for reproducibility
flu <- flu[order(flu$epitope, flu$subject, flu$cdr3b), ]
rownames(flu) <- NULL

cat("flu dataset:", nrow(flu), "paired TCRs,",
    length(unique(flu$epitope)), "epitopes\n")
print(table(flu$epitope))

usethis::use_data(flu, overwrite = TRUE, compress = "xz")
