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


# ---- as_tcr_df tests -------------------------------------------------------

test_that("as_tcr_df passes through canonical columns", {
    df <- data.frame(
        va = "TRAV1-1", cdr3a = "CAVRD", vb = "TRBV5-1", cdr3b = "CASS",
        extra = "keep_me", stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_equal(result$va, "TRAV1-1")
    expect_equal(result$cdr3b, "CASS")
    expect_true("extra" %in% colnames(result))
})

test_that("as_tcr_df detects canonical format", {
    df <- data.frame(
        va = "TRAV1-1", cdr3a = "CAVRD", vb = "TRBV5-1", cdr3b = "CASS",
        stringsAsFactors = FALSE
    )
    expect_message(as_tcr_df(df), "canonical")
})

test_that("as_tcr_df with custom col_map", {
    df <- data.frame(
        alpha_v = "TRAV1-1", alpha_cdr3 = "CAVRD",
        beta_v = "TRBV5-1", beta_cdr3 = "CASS",
        sample_id = "S1", stringsAsFactors = FALSE
    )
    result <- suppressMessages(
        as_tcr_df(df, col_map = c(va = "alpha_v", cdr3a = "alpha_cdr3",
                                   vb = "beta_v", cdr3b = "beta_cdr3"),
                  normalize_genes = FALSE)
    )
    expect_equal(result$va, "TRAV1-1")
    expect_equal(result$cdr3a, "CAVRD")
    expect_equal(result$vb, "TRBV5-1")
    expect_equal(result$cdr3b, "CASS")
    expect_true("sample_id" %in% colnames(result))
})

test_that("as_tcr_df detects scirpy/dandelion format", {
    df <- data.frame(
        IR_VJ_1_v_call = c("TRAV1-1", "TRAV1-2"),
        IR_VJ_1_j_call = c("TRAJ33", "TRAJ20"),
        IR_VJ_1_junction_aa = c("CAVRD", "CAVKD"),
        IR_VDJ_1_v_call = c("TRBV5-1", "TRBV6-1"),
        IR_VDJ_1_j_call = c("TRBJ2-7", "TRBJ1-1"),
        IR_VDJ_1_junction_aa = c("CASSIR", "CASSIK"),
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_equal(result$va, c("TRAV1-1", "TRAV1-2"))
    expect_equal(result$ja, c("TRAJ33", "TRAJ20"))
    expect_equal(result$cdr3a, c("CAVRD", "CAVKD"))
    expect_equal(result$vb, c("TRBV5-1", "TRBV6-1"))
    expect_equal(result$jb, c("TRBJ2-7", "TRBJ1-1"))
    expect_equal(result$cdr3b, c("CASSIR", "CASSIK"))
    expect_message(as_tcr_df(df), "scirpy")
})

test_that("as_tcr_df with format='scirpy' forced", {
    df <- data.frame(
        IR_VJ_1_v_call = "TRAV1-1",
        IR_VJ_1_junction_aa = "CAVRD",
        IR_VDJ_1_v_call = "TRBV5-1",
        IR_VDJ_1_junction_aa = "CASS",
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(
        as_tcr_df(df, format = "scirpy", normalize_genes = FALSE)
    )
    expect_equal(result$va, "TRAV1-1")
    expect_equal(result$cdr3b, "CASS")
})

test_that("as_tcr_df detects tcrdist3 format", {
    df <- data.frame(
        v_a_gene = "TRAV1-1", cdr3_a_aa = "CAVRD",
        v_b_gene = "TRBV5-1", cdr3_b_aa = "CASS",
        j_a_gene = "TRAJ33", j_b_gene = "TRBJ2-7",
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_equal(result$va, "TRAV1-1")
    expect_equal(result$ja, "TRAJ33")
    expect_equal(result$cdr3b, "CASS")
    expect_message(as_tcr_df(df), "tcrdist3")
})

test_that("as_tcr_df detects scRepertoire format", {
    df <- data.frame(
        CTgene = c(
            "TRAV1-1.TRAJ33.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
            "TRAV1-2.TRAJ20.TRAC_TRBV6-1.TRBD1.TRBJ1-1.TRBC1"
        ),
        CTaa = c(
            "CAVRDSSYKLIF_CASSIRSSYEQYF",
            "CAVKDSSYKLIF_CASSIKSSYEQYF"
        ),
        sample = c("S1", "S2"),
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_equal(result$va, c("TRAV1-1", "TRAV1-2"))
    expect_equal(result$ja, c("TRAJ33", "TRAJ20"))
    expect_equal(result$cdr3a, c("CAVRDSSYKLIF", "CAVKDSSYKLIF"))
    expect_equal(result$vb, c("TRBV5-1", "TRBV6-1"))
    expect_equal(result$jb, c("TRBJ2-7", "TRBJ1-1"))
    expect_equal(result$cdr3b, c("CASSIRSSYEQYF", "CASSIKSSYEQYF"))
    expect_true("sample" %in% colnames(result))
    expect_message(as_tcr_df(df), "screpertoire")
})

test_that("as_tcr_df handles scRepertoire with _TCR suffix", {
    df <- data.frame(
        CTgene_TCR = "TRAV1-1.TRAJ33.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
        CTaa_TCR = "CAVRDSSYKLIF_CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_equal(result$va, "TRAV1-1")
    expect_equal(result$cdr3b, "CASSIRSSYEQYF")
})

test_that("as_tcr_df handles scRepertoire with CTnt (nucleotide)", {
    df <- data.frame(
        CTgene = "TRAV1-1.TRAJ33.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
        CTaa = "CAVRD_CASSIR",
        CTnt = "TGTGCTGTG_TGCGCTAGC",
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_true("cdr3a_nucseq" %in% colnames(result))
    expect_true("cdr3b_nucseq" %in% colnames(result))
    expect_equal(result$cdr3a_nucseq, "TGTGCTGTG")
    expect_equal(result$cdr3b_nucseq, "TGCGCTAGC")
})

test_that("as_tcr_df handles scRepertoire multi-chain (semicolon)", {
    df <- data.frame(
        CTgene = "TRAV1-1.TRAJ33.TRAC;TRAV2.TRAJ10.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
        CTaa = "CAVRD;CAVKD_CASSIR",
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    # Should use first chain
    expect_equal(result$va, "TRAV1-1")
    expect_equal(result$cdr3a, "CAVRD")
})

test_that("as_tcr_df normalize_genes appends *01", {
    df <- data.frame(
        va = "TRAV1-1", cdr3a = "CAVRD", vb = "TRBV5-1*01", cdr3b = "CASS",
        ja = "TRAJ33", stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = TRUE))
    expect_equal(result$va, "TRAV1-1*01")
    expect_equal(result$vb, "TRBV5-1*01")
    expect_equal(result$ja, "TRAJ33*01")
})

test_that("as_tcr_df drop_incomplete removes NA rows", {
    df <- data.frame(
        va = c("TRAV1-1", NA), cdr3a = c("CAVRD", "CAVKD"),
        vb = c("TRBV5-1", "TRBV6-1"), cdr3b = c("CASS", "CASK"),
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_equal(nrow(result), 1L)
})

test_that("as_tcr_df drop_incomplete removes empty string rows", {
    df <- data.frame(
        va = c("TRAV1-1", ""), cdr3a = c("CAVRD", "CAVKD"),
        vb = c("TRBV5-1", "TRBV6-1"), cdr3b = c("CASS", "CASK"),
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_equal(nrow(result), 1L)
})

test_that("as_tcr_df drop_incomplete = FALSE preserves NAs", {
    df <- data.frame(
        va = c("TRAV1-1", NA), cdr3a = c("CAVRD", "CAVKD"),
        vb = c("TRBV5-1", "TRBV6-1"), cdr3b = c("CASS", "CASK"),
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(
        as_tcr_df(df, normalize_genes = FALSE, drop_incomplete = FALSE)
    )
    expect_equal(nrow(result), 2L)
})

test_that("as_tcr_df errors on non-data.frame input", {
    expect_error(as_tcr_df("not_a_df"), "must be a data.frame")
})

test_that("as_tcr_df errors on empty data.frame", {
    df <- data.frame(va = character(0))
    expect_error(as_tcr_df(df), "zero rows")
})

test_that("as_tcr_df errors on unrecognized columns without col_map", {
    df <- data.frame(x = "a", y = "b", stringsAsFactors = FALSE)
    expect_error(as_tcr_df(df), "auto-detect")
})

test_that("as_tcr_df errors on bad col_map format", {
    df <- data.frame(x = "a", stringsAsFactors = FALSE)
    expect_error(as_tcr_df(df, col_map = c("x", "y")), "named character")
})

test_that("as_tcr_df coerces factors to character", {
    df <- data.frame(
        va = factor("TRAV1-1"), cdr3a = factor("CAVRD"),
        vb = factor("TRBV5-1"), cdr3b = factor("CASS")
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_true(is.character(result$va))
    expect_true(is.character(result$cdr3a))
})

# ---- scRepertoire parser unit tests -----------------------------------------

test_that(".parse_screpertoire_ctgene handles standard format", {
    genes <- .parse_screpertoire_ctgene(
        "TRAV1-1.TRAJ33.TRAC_TRBV5-1.TRBD1.TRBJ2-7.TRBC2"
    )
    expect_equal(genes$va, "TRAV1-1")
    expect_equal(genes$ja, "TRAJ33")
    expect_equal(genes$vb, "TRBV5-1")
    expect_equal(genes$jb, "TRBJ2-7")
})

test_that(".parse_screpertoire_ctgene handles beta V.J.C (no D)", {
    genes <- .parse_screpertoire_ctgene(
        "TRAV1-1.TRAJ33.TRAC_TRBV5-1.TRBJ2-7.TRBC2"
    )
    expect_equal(genes$vb, "TRBV5-1")
    expect_equal(genes$jb, "TRBJ2-7")
})

test_that(".parse_screpertoire_ctgene handles NA and empty", {
    genes <- .parse_screpertoire_ctgene(c(NA, "", "NA"))
    expect_true(all(is.na(genes$va)))
    expect_true(all(is.na(genes$vb)))
})

test_that(".parse_screpertoire_cdr3 handles standard format", {
    cdr3 <- .parse_screpertoire_cdr3("CAVRDSSYKLIF_CASSIRSSYEQYF")
    expect_equal(cdr3$alpha, "CAVRDSSYKLIF")
    expect_equal(cdr3$beta, "CASSIRSSYEQYF")
})

test_that(".parse_screpertoire_cdr3 handles NA", {
    cdr3 <- .parse_screpertoire_cdr3(c(NA, ""))
    expect_true(all(is.na(cdr3$alpha)))
    expect_true(all(is.na(cdr3$beta)))
})
