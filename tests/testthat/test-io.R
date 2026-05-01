# Tests for TCR data I/O: readers, format detection, and normalization.


# ===========================================================================
# .validate_tcr_df (used by all readers)
# ===========================================================================

test_that(".validate_tcr_df detects missing columns, coerces factors, rejects NAs", {
    # Missing columns
    df <- data.frame(va = "X", cdr3a = "Y", stringsAsFactors = FALSE)
    expect_error(.validate_tcr_df(df, chains = "AB"), "missing required columns")
    expect_silent(.validate_tcr_df(df, chains = "A"))
    # Factor coercion
    df_fac <- data.frame(va = factor("TRAV1-1*01"), cdr3a = factor("CAVRD"),
                          stringsAsFactors = FALSE)
    result <- .validate_tcr_df(df_fac, chains = "A")
    expect_true(is.character(result$va))
    expect_true(is.character(result$cdr3a))
    # NA rejection
    df_na <- data.frame(va = NA_character_, cdr3a = "CAVRD",
                         stringsAsFactors = FALSE)
    expect_error(.validate_tcr_df(df_na, chains = "A"), "NA values")
})


# ===========================================================================
# read_tcr_table (generic)
# ===========================================================================

test_that("read_tcr_table reads generic TSV with auto-detection and options", {
    f <- test_path("fixtures", "test_generic.tsv")
    df <- read_tcr_table(f)
    expect_true(is.data.frame(df))
    expect_true(all(c("va", "cdr3a", "vb", "cdr3b", "ja", "jb") %in%
                        colnames(df)))
    expect_equal(nrow(df), 5L)

    # Explicit col_map
    cm <- list(va = "v_a_gene", cdr3a = "cdr3_a_aa",
               vb = "v_b_gene", cdr3b = "cdr3_b_aa")
    df2 <- read_tcr_table(f, col_map = cm)
    expect_true(all(c("va", "cdr3a", "vb", "cdr3b") %in% colnames(df2)))

    # Gene normalization
    df3 <- read_tcr_table(f, normalize_genes = TRUE)
    expect_true(all(grepl("\\*", df3$va)))
    expect_true(all(grepl("\\*", df3$vb)))

    # TCRrep compatibility
    obj <- TCRrep(df, organism = "human")
    expect_s4_class(obj, "TCRrep")

    # Error on missing file
    expect_error(read_tcr_table("nonexistent_file.tsv"), "file not found")
})


# ===========================================================================
# read_airr
# ===========================================================================

test_that("read_airr reads fixture with filtering and normalization", {
    f <- test_path("fixtures", "test_airr.tsv")
    df <- read_airr(f)
    expect_true(is.data.frame(df))
    expect_true(all(c("cell_id", "va", "ja", "cdr3a",
                       "vb", "jb", "cdr3b") %in% colnames(df)))

    # Productive filtering keeps cell4
    df_prod <- read_airr(f, productive_only = TRUE)
    expect_equal(nrow(df_prod), 5L)
    expect_true("cell4" %in% df_prod$cell_id)

    # Without filter keeps more
    df_all <- read_airr(f, productive_only = FALSE)
    expect_true(nrow(df_all) >= 5L)

    # Normalization
    df_norm <- read_airr(f, normalize_genes = TRUE)
    expect_true(all(grepl("\\*", df_norm$va)))
    expect_true(all(grepl("\\*", df_norm$vb)))

    # TCRrep compatibility
    expect_s4_class(TCRrep(df, organism = "human"), "TCRrep")

    # Error on missing file
    expect_error(read_airr("nonexistent.tsv"), "file not found")
})


# ===========================================================================
# read_10x
# ===========================================================================

test_that("read_10x reads fixture with filtering and deduplication", {
    f <- test_path("fixtures", "test_10x.csv")
    df <- read_10x(f)
    expect_true(is.data.frame(df))
    expect_true(all(c("barcode", "va", "ja", "cdr3a",
                       "vb", "jb", "cdr3b") %in% colnames(df)))

    # Filters non-cell barcodes
    expect_equal(nrow(df), 5L)
    expect_false("FFFF-1" %in% df$barcode)

    # Deduplicates by highest UMI
    dddd <- df[df$barcode == "DDDD-1", ]
    expect_equal(nrow(dddd), 1L)
    expect_equal(dddd$va, "TRAV14DV4*01")

    # Normalization
    df_norm <- read_10x(f, normalize_genes = TRUE)
    expect_true(all(grepl("\\*", df_norm$va)))

    # TCRrep compatibility
    expect_s4_class(TCRrep(df, organism = "human"), "TCRrep")

    # Error on missing file
    expect_error(read_10x("nonexistent.csv"), "file not found")
})


# ===========================================================================
# read_adaptive
# ===========================================================================

test_that("read_adaptive reads fixture with gene conversion", {
    f <- test_path("fixtures", "test_adaptive.tsv")
    df <- read_adaptive(f)
    expect_true(is.data.frame(df))
    expect_true(all(c("vb", "jb", "cdr3b") %in% colnames(df)))
    expect_equal(nrow(df), 5L)

    # Gene name conversion
    df_norm <- read_adaptive(f, normalize_genes = TRUE)
    expect_equal(df_norm$vb, c("TRBV19*01", "TRBV20-1*01", "TRBV5-1*01",
                                "TRBV28*01", "TRBV7-2*01"))
    expect_equal(df_norm$jb, c("TRBJ2-7*01", "TRBJ2-1*01", "TRBJ1-1*01",
                                "TRBJ2-3*01", "TRBJ1-5*01"))

    # Count column
    expect_true("count" %in% colnames(df))
    expect_equal(df$count, c(15L, 8L, 12L, 5L, 3L))

    # Error on missing file
    expect_error(read_adaptive("nonexistent.tsv"), "file not found")
})

test_that(".convert_adaptive_gene handles prefixes, leading zeros, and NA", {
    result <- .convert_adaptive_gene(
        c("TCRBV19*01", "TCRBJ02-07*01", "TCRAV12-02*01",
          "TCRBV05-01*01", "TCRBV07-02*01", "TCRBJ01-05*01",
          NA_character_, "", "unresolved", "TCRBV19*01")
    )
    expect_equal(result[1], "TRBV19*01")
    expect_equal(result[2], "TRBJ2-7*01")
    expect_equal(result[3], "TRAV12-2*01")
    expect_equal(result[4], "TRBV5-1*01")
    expect_equal(result[5], "TRBV7-2*01")
    expect_equal(result[6], "TRBJ1-5*01")
    expect_true(is.na(result[7]))
    expect_true(is.na(result[8]))
    expect_true(is.na(result[9]))
    expect_equal(result[10], "TRBV19*01")
})


# ===========================================================================
# as_tcr_df: format detection and conversion
# ===========================================================================

test_that("as_tcr_df handles canonical format with options", {
    df <- data.frame(
        va = "TRAV1-1", cdr3a = "CAVRD", vb = "TRBV5-1", cdr3b = "CASS",
        extra = "keep_me", stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_equal(result$va, "TRAV1-1")
    expect_equal(result$cdr3b, "CASS")
    expect_true("extra" %in% colnames(result))
    expect_message(as_tcr_df(df), "canonical")

    # Custom col_map
    df2 <- data.frame(
        alpha_v = "TRAV1-1", alpha_cdr3 = "CAVRD",
        beta_v = "TRBV5-1", beta_cdr3 = "CASS",
        sample_id = "S1", stringsAsFactors = FALSE
    )
    result2 <- suppressMessages(
        as_tcr_df(df2, col_map = c(va = "alpha_v", cdr3a = "alpha_cdr3",
                                     vb = "beta_v", cdr3b = "beta_cdr3"),
                  normalize_genes = FALSE)
    )
    expect_equal(result2$va, "TRAV1-1")
    expect_true("sample_id" %in% colnames(result2))

    # Gene normalization
    df3 <- data.frame(
        va = "TRAV1-1", cdr3a = "CAVRD", vb = "TRBV5-1*01", cdr3b = "CASS",
        ja = "TRAJ33", stringsAsFactors = FALSE
    )
    result3 <- suppressMessages(as_tcr_df(df3, normalize_genes = TRUE))
    expect_equal(result3$va, "TRAV1-1*01")
    expect_equal(result3$ja, "TRAJ33*01")

    # Factor coercion
    df_fac <- data.frame(
        va = factor("TRAV1-1"), cdr3a = factor("CAVRD"),
        vb = factor("TRBV5-1"), cdr3b = factor("CASS")
    )
    result_fac <- suppressMessages(as_tcr_df(df_fac, normalize_genes = FALSE))
    expect_true(is.character(result_fac$va))
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
    expect_equal(result$cdr3b, c("CASSIR", "CASSIK"))
    expect_message(as_tcr_df(df), "scirpy")

    # Forced format
    df2 <- data.frame(
        IR_VJ_1_v_call = "TRAV1-1", IR_VJ_1_junction_aa = "CAVRD",
        IR_VDJ_1_v_call = "TRBV5-1", IR_VDJ_1_junction_aa = "CASS",
        stringsAsFactors = FALSE
    )
    result2 <- suppressMessages(
        as_tcr_df(df2, format = "scirpy", normalize_genes = FALSE)
    )
    expect_equal(result2$va, "TRAV1-1")
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
    expect_equal(result$cdr3b, "CASS")
    expect_message(as_tcr_df(df), "tcrdist3")
})

test_that("as_tcr_df detects scRepertoire format with variants", {
    # Standard
    df <- data.frame(
        CTgene = c(
            "TRAV1-1.TRAJ33.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
            "TRAV1-2.TRAJ20.TRAC_TRBV6-1.TRBD1.TRBJ1-1.TRBC1"
        ),
        CTaa = c("CAVRDSSYKLIF_CASSIRSSYEQYF",
                  "CAVKDSSYKLIF_CASSIKSSYEQYF"),
        sample = c("S1", "S2"),
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(as_tcr_df(df, normalize_genes = FALSE))
    expect_equal(result$va, c("TRAV1-1", "TRAV1-2"))
    expect_equal(result$cdr3b, c("CASSIRSSYEQYF", "CASSIKSSYEQYF"))
    expect_true("sample" %in% colnames(result))
    expect_message(as_tcr_df(df), "screpertoire")

    # _TCR suffix variant
    df2 <- data.frame(
        CTgene_TCR = "TRAV1-1.TRAJ33.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
        CTaa_TCR = "CAVRDSSYKLIF_CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    result2 <- suppressMessages(as_tcr_df(df2, normalize_genes = FALSE))
    expect_equal(result2$va, "TRAV1-1")

    # CTnt nucleotide variant
    df3 <- data.frame(
        CTgene = "TRAV1-1.TRAJ33.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
        CTaa = "CAVRD_CASSIR",
        CTnt = "TGTGCTGTG_TGCGCTAGC",
        stringsAsFactors = FALSE
    )
    result3 <- suppressMessages(as_tcr_df(df3, normalize_genes = FALSE))
    expect_equal(result3$cdr3a_nucseq, "TGTGCTGTG")
    expect_equal(result3$cdr3b_nucseq, "TGCGCTAGC")

    # Multi-chain semicolon variant
    df4 <- data.frame(
        CTgene = "TRAV1-1.TRAJ33.TRAC;TRAV2.TRAJ10.TRAC_TRBV5-1.None.TRBJ2-7.TRBC2",
        CTaa = "CAVRD;CAVKD_CASSIR",
        stringsAsFactors = FALSE
    )
    result4 <- suppressMessages(as_tcr_df(df4, normalize_genes = FALSE))
    expect_equal(result4$va, "TRAV1-1")
    expect_equal(result4$cdr3a, "CAVRD")
})

test_that("as_tcr_df drop_incomplete and error handling", {
    # NA rows dropped by default
    df_na <- data.frame(
        va = c("TRAV1-1", NA), cdr3a = c("CAVRD", "CAVKD"),
        vb = c("TRBV5-1", "TRBV6-1"), cdr3b = c("CASS", "CASK"),
        stringsAsFactors = FALSE
    )
    expect_equal(nrow(suppressMessages(as_tcr_df(df_na))), 1L)

    # Empty string rows dropped
    df_empty <- data.frame(
        va = c("TRAV1-1", ""), cdr3a = c("CAVRD", "CAVKD"),
        vb = c("TRBV5-1", "TRBV6-1"), cdr3b = c("CASS", "CASK"),
        stringsAsFactors = FALSE
    )
    expect_equal(nrow(suppressMessages(as_tcr_df(df_empty))), 1L)

    # drop_incomplete=FALSE preserves NAs
    result <- suppressMessages(
        as_tcr_df(df_na, normalize_genes = FALSE, drop_incomplete = FALSE)
    )
    expect_equal(nrow(result), 2L)

    # Error cases
    expect_error(as_tcr_df("not_a_df"), "must be a data.frame")
    expect_error(as_tcr_df(data.frame(va = character(0))), "zero rows")
    expect_error(as_tcr_df(data.frame(x = "a", y = "b",
                                       stringsAsFactors = FALSE)),
                 "auto-detect")
    expect_error(as_tcr_df(data.frame(x = "a", stringsAsFactors = FALSE),
                            col_map = c("x", "y")),
                 "named character")
})


# ===========================================================================
# scRepertoire parser internals
# ===========================================================================

test_that(".parse_screpertoire_ctgene handles standard, no-D, and NA", {
    genes <- .parse_screpertoire_ctgene(
        "TRAV1-1.TRAJ33.TRAC_TRBV5-1.TRBD1.TRBJ2-7.TRBC2"
    )
    expect_equal(genes$va, "TRAV1-1")
    expect_equal(genes$ja, "TRAJ33")
    expect_equal(genes$vb, "TRBV5-1")
    expect_equal(genes$jb, "TRBJ2-7")

    # No D gene variant
    genes2 <- .parse_screpertoire_ctgene(
        "TRAV1-1.TRAJ33.TRAC_TRBV5-1.TRBJ2-7.TRBC2"
    )
    expect_equal(genes2$vb, "TRBV5-1")
    expect_equal(genes2$jb, "TRBJ2-7")

    # NA/empty
    genes3 <- .parse_screpertoire_ctgene(c(NA, "", "NA"))
    expect_true(all(is.na(genes3$va)))
    expect_true(all(is.na(genes3$vb)))
})

test_that(".parse_screpertoire_cdr3 handles standard and NA", {
    cdr3 <- .parse_screpertoire_cdr3("CAVRDSSYKLIF_CASSIRSSYEQYF")
    expect_equal(cdr3$alpha, "CAVRDSSYKLIF")
    expect_equal(cdr3$beta, "CASSIRSSYEQYF")

    cdr3_na <- .parse_screpertoire_cdr3(c(NA, ""))
    expect_true(all(is.na(cdr3_na$alpha)))
    expect_true(all(is.na(cdr3_na$beta)))
})
