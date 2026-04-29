# ---- io_utils.R tests -------------------------------------------------------

test_that(".normalize_gene_names appends *01 when missing", {
    input    <- c("TRAV1-1", "TRBV19*01", "TRAV12-2*02", "TRBJ2-7")
    expected <- c("TRAV1-1*01", "TRBV19*01", "TRAV12-2*02", "TRBJ2-7*01")
    expect_equal(.normalize_gene_names(input), expected)
})

test_that(".normalize_gene_names handles NA and empty strings", {
    input    <- c(NA_character_, "", "TRAV1-1")
    expected <- c(NA_character_, "", "TRAV1-1*01")
    expect_equal(.normalize_gene_names(input), expected)
})

test_that(".validate_tcr_df detects missing columns", {
    df <- data.frame(va = "X", cdr3a = "Y", stringsAsFactors = FALSE)
    expect_error(.validate_tcr_df(df, chains = "AB"), "missing required columns")
    expect_silent(.validate_tcr_df(df, chains = "A"))
})

test_that(".validate_tcr_df coerces factors", {
    df <- data.frame(va = factor("TRAV1-1*01"), cdr3a = factor("CAVRD"),
                     stringsAsFactors = FALSE)
    result <- .validate_tcr_df(df, chains = "A")
    expect_true(is.character(result$va))
    expect_true(is.character(result$cdr3a))
})

test_that(".validate_tcr_df rejects NAs in required columns", {
    df <- data.frame(va = NA_character_, cdr3a = "CAVRD",
                     stringsAsFactors = FALSE)
    expect_error(.validate_tcr_df(df, chains = "A"), "NA values")
})

test_that(".standardize_columns renames matching columns", {
    df <- data.frame(v_a_gene = "X", cdr3_a_aa = "Y", other = "Z",
                     stringsAsFactors = FALSE)
    col_map <- list(va = "v_a_gene", cdr3a = "cdr3_a_aa")
    result <- .standardize_columns(df, col_map)
    expect_true("va" %in% colnames(result))
    expect_true("cdr3a" %in% colnames(result))
    expect_true("other" %in% colnames(result))
})

test_that(".standardize_columns skips missing source columns", {
    df <- data.frame(va = "X", stringsAsFactors = FALSE)
    col_map <- list(va = "v_a_gene", cdr3a = "cdr3_a_aa")
    result <- .standardize_columns(df, col_map)
    expect_equal(colnames(result), "va")
})

# ---- io_generic.R tests ----------------------------------------------------

test_that("read_tcr_table reads tcrdist3-style TSV with auto-detection", {
    f <- test_path("fixtures", "test_generic.tsv")
    df <- read_tcr_table(f)
    expect_true(is.data.frame(df))
    expect_true(all(c("va", "cdr3a", "vb", "cdr3b") %in% colnames(df)))
    expect_equal(nrow(df), 5L)
})

test_that("read_tcr_table renames ja/jb columns", {
    f <- test_path("fixtures", "test_generic.tsv")
    df <- read_tcr_table(f)
    expect_true("ja" %in% colnames(df))
    expect_true("jb" %in% colnames(df))
})

test_that("read_tcr_table with explicit col_map", {
    f <- test_path("fixtures", "test_generic.tsv")
    cm <- list(va = "v_a_gene", cdr3a = "cdr3_a_aa",
               vb = "v_b_gene", cdr3b = "cdr3_b_aa")
    df <- read_tcr_table(f, col_map = cm)
    expect_true(all(c("va", "cdr3a", "vb", "cdr3b") %in% colnames(df)))
    expect_equal(nrow(df), 5L)
})

test_that("read_tcr_table errors on missing file", {
    expect_error(read_tcr_table("nonexistent_file.tsv"), "file not found")
})

test_that("read_tcr_table gene names are normalized", {
    f <- test_path("fixtures", "test_generic.tsv")
    df <- read_tcr_table(f, normalize_genes = TRUE)
    # All gene names should already contain *
    expect_true(all(grepl("\\*", df$va)))
    expect_true(all(grepl("\\*", df$vb)))
})

test_that("read_tcr_table output compatible with TCRrep", {
    f <- test_path("fixtures", "test_generic.tsv")
    df <- read_tcr_table(f)
    obj <- TCRrep(df, organism = "human")
    expect_s4_class(obj, "TCRrep")
    expect_equal(nrow(obj@clone_df), 5L)
})

# ---- io_airr.R tests -------------------------------------------------------

test_that("read_airr reads AIRR fixture and returns paired data", {
    f <- test_path("fixtures", "test_airr.tsv")
    df <- read_airr(f)
    expect_true(is.data.frame(df))
    expect_true(all(c("cell_id", "va", "ja", "cdr3a",
                       "vb", "jb", "cdr3b") %in% colnames(df)))
})

test_that("read_airr filters non-productive and still pairs cell4", {
    f <- test_path("fixtures", "test_airr.tsv")
    df <- read_airr(f, productive_only = TRUE)
    # cell4 has one non-productive TRA (seq7) filtered, one productive (seq8) kept
    expect_equal(nrow(df), 5L)
    expect_true("cell4" %in% df$cell_id)
})

test_that("read_airr without productive filter keeps all", {
    f <- test_path("fixtures", "test_airr.tsv")
    df <- read_airr(f, productive_only = FALSE)
    # cell4 has 2 TRA rows when not filtering -> merge duplicates
    expect_true(nrow(df) >= 5L)
})

test_that("read_airr normalizes gene names", {
    f <- test_path("fixtures", "test_airr.tsv")
    df <- read_airr(f, normalize_genes = TRUE)
    # All gene names should contain * allele suffix
    expect_true(all(grepl("\\*", df$va)))
    expect_true(all(grepl("\\*", df$vb)))
    expect_true(all(grepl("\\*", df$ja)))
    expect_true(all(grepl("\\*", df$jb)))
})

test_that("read_airr errors on missing file", {
    expect_error(read_airr("nonexistent.tsv"), "file not found")
})

test_that("read_airr output compatible with TCRrep", {
    f <- test_path("fixtures", "test_airr.tsv")
    df <- read_airr(f)
    obj <- TCRrep(df, organism = "human")
    expect_s4_class(obj, "TCRrep")
})

# ---- io_10x.R tests --------------------------------------------------------

test_that("read_10x reads fixture and returns paired data", {
    f <- test_path("fixtures", "test_10x.csv")
    df <- read_10x(f)
    expect_true(is.data.frame(df))
    expect_true(all(c("barcode", "va", "ja", "cdr3a",
                       "vb", "jb", "cdr3b") %in% colnames(df)))
})

test_that("read_10x filters out non-cell barcodes", {
    f <- test_path("fixtures", "test_10x.csv")
    df <- read_10x(f)
    # FFFF-1 has is_cell=False -> excluded; 5 valid barcodes remain
    expect_equal(nrow(df), 5L)
    expect_false("FFFF-1" %in% df$barcode)
})

test_that("read_10x deduplicates multi-chain by highest UMI", {
    f <- test_path("fixtures", "test_10x.csv")
    df <- read_10x(f)
    # DDDD-1 has 2 TRA contigs: TRAV14DV4*01 (umis=4) and TRAV38-2DV8*01 (umis=2)
    # Should keep higher UMI -> TRAV14DV4*01
    dddd <- df[df$barcode == "DDDD-1", ]
    expect_equal(nrow(dddd), 1L)
    expect_equal(dddd$va, "TRAV14DV4*01")
})

test_that("read_10x normalizes gene names", {
    f <- test_path("fixtures", "test_10x.csv")
    df <- read_10x(f, normalize_genes = TRUE)
    expect_true(all(grepl("\\*", df$va)))
    expect_true(all(grepl("\\*", df$vb)))
})

test_that("read_10x errors on missing file", {
    expect_error(read_10x("nonexistent.csv"), "file not found")
})

test_that("read_10x output compatible with TCRrep", {
    f <- test_path("fixtures", "test_10x.csv")
    df <- read_10x(f)
    obj <- TCRrep(df, organism = "human")
    expect_s4_class(obj, "TCRrep")
})

# ---- io_adaptive.R tests ---------------------------------------------------

test_that(".convert_adaptive_gene converts prefixes correctly", {
    input <- c("TCRBV19*01", "TCRBJ02-07*01", "TCRAV12-02*01")
    result <- .convert_adaptive_gene(input)
    expect_equal(result[1], "TRBV19*01")
    expect_equal(result[2], "TRBJ2-7*01")
    expect_equal(result[3], "TRAV12-2*01")
})

test_that(".convert_adaptive_gene strips leading zeros from family and subfamily", {
    input <- c("TCRBV05-01*01", "TCRBV07-02*01", "TCRBJ01-05*01")
    result <- .convert_adaptive_gene(input)
    expect_equal(result[1], "TRBV5-1*01")
    expect_equal(result[2], "TRBV7-2*01")
    expect_equal(result[3], "TRBJ1-5*01")
})

test_that(".convert_adaptive_gene handles NA, empty, unresolved", {
    input <- c(NA_character_, "", "unresolved", "TCRBV19*01")
    result <- .convert_adaptive_gene(input)
    expect_true(is.na(result[1]))
    expect_true(is.na(result[2]))
    expect_true(is.na(result[3]))
    expect_equal(result[4], "TRBV19*01")
})

test_that("read_adaptive reads fixture correctly", {
    f <- test_path("fixtures", "test_adaptive.tsv")
    df <- read_adaptive(f)
    expect_true(is.data.frame(df))
    expect_true(all(c("vb", "jb", "cdr3b") %in% colnames(df)))
    expect_equal(nrow(df), 5L)
})

test_that("read_adaptive converts gene names to IMGT format", {
    f <- test_path("fixtures", "test_adaptive.tsv")
    df <- read_adaptive(f, normalize_genes = TRUE)
    # Verify specific conversions
    expect_equal(df$vb[1], "TRBV19*01")
    expect_equal(df$vb[2], "TRBV20-1*01")
    expect_equal(df$vb[3], "TRBV5-1*01")
    expect_equal(df$vb[4], "TRBV28*01")
    expect_equal(df$vb[5], "TRBV7-2*01")
})

test_that("read_adaptive converts J gene names correctly", {
    f <- test_path("fixtures", "test_adaptive.tsv")
    df <- read_adaptive(f, normalize_genes = TRUE)
    expect_equal(df$jb[1], "TRBJ2-7*01")
    expect_equal(df$jb[2], "TRBJ2-1*01")
    expect_equal(df$jb[3], "TRBJ1-1*01")
    expect_equal(df$jb[4], "TRBJ2-3*01")
    expect_equal(df$jb[5], "TRBJ1-5*01")
})

test_that("read_adaptive includes count column", {
    f <- test_path("fixtures", "test_adaptive.tsv")
    df <- read_adaptive(f)
    expect_true("count" %in% colnames(df))
    expect_equal(df$count, c(15L, 8L, 12L, 5L, 3L))
})

test_that("read_adaptive errors on missing file", {
    expect_error(read_adaptive("nonexistent.tsv"), "file not found")
})
