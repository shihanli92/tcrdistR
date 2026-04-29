test_that("TCRrep basic creation succeeds and returns valid object", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human")
    expect_true(is(obj, "TCRrep"))
    expect_equal(nrow(obj@clone_df), 3L)
    expect_equal(obj@organism, "human")
    expect_equal(obj@chains,   "AB")
    expect_equal(obj@metric,   "tcrdist")
})

test_that("TCRrep rejects invalid chains argument", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF"),
        vb    = c("TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )

    expect_error(TCRrep(tcrs, "human", chains = "XYZ"))
})

test_that("TCRrep rejects clone_df with missing required columns for chains='AB'", {
    # Only alpha columns — missing vb and cdr3b
    tcrs_alpha_only <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        stringsAsFactors = FALSE
    )

    expect_error(TCRrep(tcrs_alpha_only, "human", chains = "AB"))
})

test_that("TCRrep single-chain 'A' works with only alpha columns", {
    tcrs_alpha <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs_alpha, "human", chains = "A")
    expect_true(is(obj, "TCRrep"))
    expect_equal(obj@chains, "A")
    expect_equal(nrow(obj@clone_df), 3L)
})

test_that("TCRrep coerces factor columns to character", {
    tcrs_factor <- data.frame(
        va    = factor(c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01")),
        cdr3a = factor(c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF")),
        vb    = factor(c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01")),
        cdr3b = factor(c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"))
    )

    obj <- TCRrep(tcrs_factor, "human")
    expect_true(is(obj, "TCRrep"))
    expect_true(is.character(obj@clone_df$va))
    expect_true(is.character(obj@clone_df$cdr3a))
    expect_true(is.character(obj@clone_df$vb))
    expect_true(is.character(obj@clone_df$cdr3b))
})

test_that("TCRrep with compute_distances=TRUE stores a matrix in paired_dist", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", compute_distances = TRUE)
    expect_false(is.null(obj@paired_dist))
    expect_true(is.matrix(obj@paired_dist))
    expect_equal(dim(obj@paired_dist), c(3L, 3L))
})

test_that("TCRrep with compute_distances=FALSE (default) leaves paired_dist NULL", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human")
    expect_true(is.null(obj@paired_dist))
})

test_that("show method produces expected output for TCRrep", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human")
    out <- capture.output(show(obj))
    expect_true(any(grepl("TCRrep object", out)))
    expect_true(any(grepl("clonotypes", out)))
})

test_that("TCRrep succeeds with empty data.frame", {
    empty_df <- data.frame(
        va    = character(0L),
        cdr3a = character(0L),
        vb    = character(0L),
        cdr3b = character(0L),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(empty_df, "human")
    expect_true(is(obj, "TCRrep"))
    expect_equal(nrow(obj@clone_df), 0L)
})

test_that("TCRrep default weights and gap penalties match package constants", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF"),
        vb    = c("TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human")
    expect_equal(obj@weights$cdr3,          3L)
    expect_equal(obj@weights$v_region,       1L)
    expect_equal(obj@gap_penalties$cdr3,    12L)
    expect_equal(obj@gap_penalties$v_region,  4L)
})
