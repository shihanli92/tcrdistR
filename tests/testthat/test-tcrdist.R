test_that("tcrdist_matrix returns correct dimensions", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
        stringsAsFactors = FALSE
    )
    mat <- tcrdist_matrix(tcrs, "human")
    expect_equal(nrow(mat), 2L)
    expect_equal(ncol(mat), 2L)
})

test_that("tcrdist_matrix is symmetric with zero diagonal", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-2*01", "TRAV12-1*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CALSDRSYEKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV5-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CASSLGTEAFF"),
        stringsAsFactors = FALSE
    )
    mat <- tcrdist_matrix(tcrs, "human")
    expect_true(isSymmetric(mat))
    expect_equal(unname(diag(mat)), rep(0.0, 3))
})

test_that("identical TCRs have distance 0", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )
    mat <- tcrdist_matrix(tcrs, "human")
    expect_equal(mat[1, 2], 0.0)
})

test_that("tcrdist_matrix handles single TCR", {
    tcrs <- data.frame(
        va    = "TRAV1-1*01",
        cdr3a = "CAVRDSSYKLIF",
        vb    = "TRBV19*01",
        cdr3b = "CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    mat <- tcrdist_matrix(tcrs, "human")
    expect_equal(nrow(mat), 1L)
    expect_equal(ncol(mat), 1L)
    expect_equal(mat[1, 1], 0.0)
})

test_that("tcrdist_matrix handles empty input", {
    tcrs <- data.frame(
        va    = character(0),
        cdr3a = character(0),
        vb    = character(0),
        cdr3b = character(0),
        stringsAsFactors = FALSE
    )
    mat <- tcrdist_matrix(tcrs, "human")
    expect_equal(nrow(mat), 0L)
    expect_equal(ncol(mat), 0L)
})

test_that("tcrdist_matrix rejects missing columns", {
    tcrs <- data.frame(va = "TRAV1-1*01", cdr3a = "CAVRDSSYKLIF",
                       stringsAsFactors = FALSE)
    expect_error(tcrdist_matrix(tcrs, "human"), "missing required columns")
})

test_that("tcrdist_matrix rejects NA values", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", NA),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )
    expect_error(tcrdist_matrix(tcrs, "human"), "NA values")
})

test_that("tcrdist_matrix rejects unknown V genes", {
    tcrs <- data.frame(
        va    = "TRAV999*01",
        cdr3a = "CAVRDSSYKLIF",
        vb    = "TRBV19*01",
        cdr3b = "CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    expect_error(tcrdist_matrix(tcrs, "human"), "not found")
})

test_that("tcrdist_matrix produces known reference values", {
    # Cross-validated against rconga: these 3 TCRs produce this exact matrix
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-2*01", "TRAV12-1*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CALSDRSYEKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV5-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CASSLGTEAFF"),
        stringsAsFactors = FALSE
    )
    mat <- tcrdist_matrix(tcrs, "human")
    expect_equal(mat[1, 2], 28)
    expect_equal(mat[1, 3], 239)
    expect_equal(mat[2, 3], 227)
})
