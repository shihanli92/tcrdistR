# Tests for R/tcr_db_matching.R — TCR database matching.


# ---------------------------------------------------------------------------
# Test 1: Column normalization
# ---------------------------------------------------------------------------

test_that(".normalize_tcr_db_columns renames va_gene to va", {
    df <- data.frame(va_gene = "X", vb_gene = "Y", cdr3a = "A",
                     stringsAsFactors = FALSE)
    result <- tcrdistR:::.normalize_tcr_db_columns(df)
    expect_true("va" %in% colnames(result))
    expect_true("vb" %in% colnames(result))
    expect_false("va_gene" %in% colnames(result))
    expect_false("vb_gene" %in% colnames(result))
})


test_that(".normalize_tcr_db_columns is a no-op if already correct", {
    df <- data.frame(va = "X", vb = "Y", stringsAsFactors = FALSE)
    result <- tcrdistR:::.normalize_tcr_db_columns(df)
    expect_identical(colnames(result), colnames(df))
})


# ---------------------------------------------------------------------------
# Test 2: Empty result scaffold
# ---------------------------------------------------------------------------

test_that(".empty_tcr_db_match_df has correct structure", {
    result <- tcrdistR:::.empty_tcr_db_match_df()
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 0L)
    expected_cols <- c("tcrdist", "pvalue_adj", "fdr_value",
                       "query_index", "db_index", "va", "cdr3a",
                       "vb", "cdr3b")
    expect_true(all(expected_cols %in% colnames(result)))
})


# ---------------------------------------------------------------------------
# Test 3: strict_single_chain_match_tcrs_to_db
# ---------------------------------------------------------------------------

test_that("strict_single_chain_match finds exact CDR3 matches in human DB", {
    skip_on_cran()

    # Use CDR3 sequences known to be in the human database
    db_path <- system.file("extdata", "human_tcr_db_for_matching.tsv",
                            package = "tcrdistR")
    skip_if(!nzchar(db_path), "human_tcr_db_for_matching.tsv not found")

    # Read a few rows from the DB to get known CDR3 sequences
    db_sample <- utils::read.delim(db_path, stringsAsFactors = FALSE,
                                    nrows = 5)
    db_sample <- tcrdistR:::.normalize_tcr_db_columns(db_sample)

    # Use the first entry's CDR3s as our query
    query_df <- data.frame(
        cdr3a = db_sample$cdr3a[1L],
        cdr3b = db_sample$cdr3b[1L],
        stringsAsFactors = FALSE
    )

    result <- strict_single_chain_match_tcrs_to_db(query_df, "human")

    expect_type(result, "list")
    expect_true("alpha_matches" %in% names(result))
    expect_true("beta_matches" %in% names(result))
    expect_s3_class(result$alpha_matches, "data.frame")
    expect_s3_class(result$beta_matches, "data.frame")

    # We should find at least 1 match for both chains since we used
    # CDR3 sequences from the database itself
    expect_true(nrow(result$alpha_matches) >= 1L)
    expect_true(nrow(result$beta_matches) >= 1L)
})


test_that("strict_single_chain_match works with mouse DB", {
    skip_on_cran()

    db_path <- system.file("extdata", "mouse_tcr_db_for_matching.tsv",
                            package = "tcrdistR")
    skip_if(!nzchar(db_path), "mouse_tcr_db_for_matching.tsv not found")

    db_sample <- utils::read.delim(db_path, stringsAsFactors = FALSE,
                                    nrows = 20)
    db_sample <- tcrdistR:::.normalize_tcr_db_columns(db_sample)

    # Mouse DB is single-chain focused: many rows have empty cdr3a.
    # Use a known non-empty cdr3b for beta chain matching.
    query_df <- data.frame(
        cdr3a = "PLACEHOLDER",
        cdr3b = db_sample$cdr3b[1L],
        stringsAsFactors = FALSE
    )

    result <- strict_single_chain_match_tcrs_to_db(query_df, "mouse")
    expect_type(result, "list")
    expect_true(nrow(result$beta_matches) >= 1L)
})


test_that("strict_single_chain_match warns for unsupported organism", {
    query_df <- data.frame(cdr3a = "CAVRD", cdr3b = "CASSI",
                           stringsAsFactors = FALSE)
    expect_warning(
        strict_single_chain_match_tcrs_to_db(query_df, "zebrafish"),
        "no default database"
    )
})


# ---------------------------------------------------------------------------
# Test 4: find_significant_tcrdist_matches input validation
# ---------------------------------------------------------------------------

test_that("find_significant_tcrdist_matches validates columns", {
    bad_query <- data.frame(va = "X", stringsAsFactors = FALSE)
    bad_db <- data.frame(va = "X", cdr3a = "Y", vb = "Z", cdr3b = "W",
                         stringsAsFactors = FALSE)

    expect_error(
        find_significant_tcrdist_matches(bad_query, bad_db, "human"),
        "missing column"
    )
})


test_that("find_significant_tcrdist_matches validates background columns", {
    query <- data.frame(
        va = "TRAV1-2*01", cdr3a = "CAVRD",
        vb = "TRBV6-4*01", cdr3b = "CASSI",
        stringsAsFactors = FALSE
    )
    db <- data.frame(
        va = "TRAV1-2*01", cdr3a = "CAVRD",
        vb = "TRBV6-4*01", cdr3b = "CASSI",
        stringsAsFactors = FALSE
    )
    # background_tcrs_df defaults to query, which is missing ja/jb/nucseqs
    expect_error(
        find_significant_tcrdist_matches(query, db, "human"),
        "missing columns"
    )
})


# ---------------------------------------------------------------------------
# Test 5: match_tcrs_to_db warns for non-human
# ---------------------------------------------------------------------------

test_that("match_tcrs_to_db warns for non-human organism without custom DB", {
    tcr_df <- data.frame(
        va = "TRAV1*01", ja = "TRAJ1*01", cdr3a = "X", cdr3a_nucseq = "x",
        vb = "TRBV1*01", jb = "TRBJ1*01", cdr3b = "Y", cdr3b_nucseq = "y",
        stringsAsFactors = FALSE
    )
    expect_warning(
        match_tcrs_to_db(tcr_df, "mouse"),
        "default paired database"
    )
})


# ---------------------------------------------------------------------------
# Test 6: find_significant_tcrdist_matches end-to-end
# ---------------------------------------------------------------------------

test_that("find_significant_tcrdist_matches end-to-end returns matches", {
    skip_on_cran()

    data(dash, envir = environment())
    # Use first 5 mouse TCRs as query (has nucseq columns for background)
    query <- dash[1:5, ]

    # Put identical TCRs in db — guarantees distance = 0 matches
    db <- data.frame(
        va    = query$va,
        cdr3a = query$cdr3a,
        vb    = query$vb,
        cdr3b = query$cdr3b,
        stringsAsFactors = FALSE
    )

    result <- find_significant_tcrdist_matches(
        query, db, "mouse",
        adjusted_pvalue_threshold = 1.0,
        num_random_samples = 500L
    )

    expect_true(is.data.frame(result))
    expect_true(nrow(result) >= 5L)
    expect_true("tcrdist" %in% colnames(result))
    expect_true("pvalue_adj" %in% colnames(result))
    expect_true("fdr_value" %in% colnames(result))
    expect_true("db_va" %in% colnames(result))
    # Self-matches should have distance 0
    self_matches <- result[result$tcrdist == 0, ]
    expect_true(nrow(self_matches) >= 5L)
    # Results sorted by pvalue_adj
    expect_equal(result$pvalue_adj, sort(result$pvalue_adj))
})

test_that("find_significant_tcrdist_matches returns empty for stringent threshold", {
    skip_on_cran()

    data(dash, envir = environment())
    query <- dash[1:3, ]

    # db with very different CDR3s (use mouse V genes, valid AAs)
    db <- data.frame(
        va    = rep("TRAV7-3*01", 3),
        cdr3a = c("CAAAAAAAAF", "CGGGGGGGGF", "CLLLLLLLLLF"),
        vb    = rep("TRBV13-1*01", 3),
        cdr3b = c("CASSSSSSSSF", "CATTTTTTTTF", "CASRRRRRRRRF"),
        stringsAsFactors = FALSE
    )

    result <- find_significant_tcrdist_matches(
        query, db, "mouse",
        adjusted_pvalue_threshold = 0.05,
        num_random_samples = 500L
    )

    expect_true(is.data.frame(result))
    expect_true(all(c("tcrdist", "pvalue_adj") %in% colnames(result)))
})
