# Tests for R/cd8_scoring.R — CD8 logistic regression scoring.


# ---------------------------------------------------------------------------
# Test 1: Model loading
# ---------------------------------------------------------------------------

test_that(".load_cd8_logreg_models loads both chains with correct structure", {
    models <- tcrdistR:::.load_cd8_logreg_models()

    expect_type(models, "list")
    expect_true("A" %in% names(models))
    expect_true("B" %in% names(models))

    for (ab in c("A", "B")) {
        m <- models[[ab]]
        expect_true("weights" %in% names(m))
        expect_true("vgene_indexer" %in% names(m))
        expect_true("jgene_indexer" %in% names(m))
        expect_true("window_size" %in% names(m))
        expect_true("min_lenbin" %in% names(m))
        expect_true("max_lenbin" %in% names(m))
        expect_true("NV" %in% names(m))
        expect_true("NJ" %in% names(m))
        expect_true("NL" %in% names(m))
        expect_true("NC" %in% names(m))

        # Weights vector length = NV + NJ + NL + NC + 1 (bias)
        NTOT <- m[["NV"]] + m[["NJ"]] + m[["NL"]] + m[["NC"]]
        expect_equal(length(m[["weights"]]), NTOT + 1L)

        # Gene indexers are 0-based
        expect_true(all(m[["vgene_indexer"]] >= 0L))
        expect_true(all(m[["jgene_indexer"]] >= 0L))

        # Model params are positive integers
        expect_true(m[["window_size"]] > 0L)
        expect_true(m[["min_lenbin"]] > 0L)
        expect_true(m[["max_lenbin"]] > m[["min_lenbin"]])
    }
})


# ---------------------------------------------------------------------------
# Test 2: Feature encoding
# ---------------------------------------------------------------------------

test_that(".encode_single_chain_tcr has correct structure", {
    models <- tcrdistR:::.load_cd8_logreg_models()
    m <- models[["A"]]

    x <- tcrdistR:::.encode_single_chain_tcr(
        "TRAV1-2*01", "TRAJ33*01", "CAVMDSSYKLIF", m
    )

    NTOT <- m[["NV"]] + m[["NJ"]] + m[["NL"]] + m[["NC"]]

    # Correct length
    expect_equal(length(x), NTOT + 1L)

    # Bias term is 1
    expect_equal(x[NTOT + 1L], 1.0)

    # Exactly one V-gene slot is 1
    NV <- length(m[["vgene_indexer"]]) + 1L
    expect_equal(sum(x[seq_len(NV)]), 1.0)

    # Exactly one J-gene slot is 1
    NJ <- length(m[["jgene_indexer"]]) + 1L
    expect_equal(sum(x[seq(NV + 1L, NV + NJ)]), 1.0)

    # Exactly one length-bin slot is 1
    NL <- m[["max_lenbin"]] - m[["min_lenbin"]] + 1L
    expect_equal(sum(x[seq(NV + NJ + 1L, NV + NJ + NL)]), 1.0)
})


test_that(".encode_single_chain_tcr uses UNK slot for unknown genes", {
    models <- tcrdistR:::.load_cd8_logreg_models()
    m <- models[["A"]]

    x <- tcrdistR:::.encode_single_chain_tcr(
        "TRAV_FAKE*01", "TRAJ_FAKE*01", "CAVMDSSYKLIF", m
    )

    NV <- length(m[["vgene_indexer"]]) + 1L
    NJ <- length(m[["jgene_indexer"]]) + 1L

    # UNK V-gene slot (last one) should be 1
    expect_equal(x[NV], 1.0)

    # UNK J-gene slot (last one) should be 1
    expect_equal(x[NV + NJ], 1.0)
})


# ---------------------------------------------------------------------------
# Test 3: Scoring function
# ---------------------------------------------------------------------------

test_that("make_cd8_score_table_column returns finite scores", {
    tcr_df <- data.frame(
        va = c("TRAV1-2*01", "TRAV12-2*01"),
        ja = c("TRAJ33*01", "TRAJ49*01"),
        cdr3a = c("CAVMDSSYKLIF", "CAVSANSGTYF"),
        vb = c("TRBV6-4*01", "TRBV20-1*01"),
        jb = c("TRBJ2-1*01", "TRBJ1-1*01"),
        cdr3b = c("CASSLAPGATNEKLFF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    scores <- make_cd8_score_table_column(tcr_df)

    expect_type(scores, "double")
    expect_equal(length(scores), 2L)
    expect_true(all(is.finite(scores)))
})


test_that("make_cd8_score_table_column with sigmoid returns values in (0,1)", {
    tcr_df <- data.frame(
        va = "TRAV1-2*01", ja = "TRAJ33*01",
        cdr3a = "CAVMDSSYKLIF",
        vb = "TRBV6-4*01", jb = "TRBJ2-1*01",
        cdr3b = "CASSLAPGATNEKLFF",
        stringsAsFactors = FALSE
    )

    scores_raw <- make_cd8_score_table_column(tcr_df, use_sigmoid = FALSE)
    scores_sig <- make_cd8_score_table_column(tcr_df, use_sigmoid = TRUE)

    expect_true(all(is.finite(scores_raw)))
    expect_true(all(scores_sig > 0 & scores_sig < 1))
})


test_that("make_cd8_score_table_column is deterministic", {
    tcr_df <- data.frame(
        va = c("TRAV1-2*01", "TRAV12-2*01"),
        ja = c("TRAJ33*01", "TRAJ49*01"),
        cdr3a = c("CAVMDSSYKLIF", "CAVSANSGTYF"),
        vb = c("TRBV6-4*01", "TRBV20-1*01"),
        jb = c("TRBJ2-1*01", "TRBJ1-1*01"),
        cdr3b = c("CASSLAPGATNEKLFF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    scores1 <- make_cd8_score_table_column(tcr_df)
    scores2 <- make_cd8_score_table_column(tcr_df)

    expect_identical(scores1, scores2)
})


# ---------------------------------------------------------------------------
# Test 4: Input validation
# ---------------------------------------------------------------------------

test_that("make_cd8_score_table_column validates inputs", {
    expect_error(
        make_cd8_score_table_column("not a data.frame"),
        "data.frame"
    )

    bad_df <- data.frame(va = "x", cdr3a = "y", stringsAsFactors = FALSE)
    expect_error(
        make_cd8_score_table_column(bad_df),
        "missing required columns"
    )
})
