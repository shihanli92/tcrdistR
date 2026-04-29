test_that("tcrdist_rect(A, A) equals tcrdist_matrix(A)", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    rect <- tcrdist_rect(tcrs, tcrs, "human")
    mat  <- tcrdist_matrix(tcrs, "human")
    expect_equal(rect, mat)
})

test_that("tcrdist_rect returns correct dimensions", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    query <- tcrs[1:2, ]
    ref   <- tcrs[2:3, ]
    rect  <- tcrdist_rect(query, ref, "human")
    expect_equal(nrow(rect), 2L)
    expect_equal(ncol(rect), 2L)
})

test_that("tcrdist_rect single query returns 1 x N matrix", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    query <- tcrs[1L, ]
    rect  <- tcrdist_rect(query, tcrs, "human")
    expect_equal(nrow(rect), 1L)
    expect_equal(ncol(rect), 3L)
})

test_that("tcrdist_rect produces known reference values", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    rect <- tcrdist_rect(tcrs, tcrs, "human")
    # Values from tcrdist_matrix on the same data
    mat <- tcrdist_matrix(tcrs, "human")
    expect_equal(rect[1, 2], mat[1, 2])
    expect_equal(rect[1, 3], mat[1, 3])
    expect_equal(rect[2, 3], mat[2, 3])
})

test_that("tcrdist_rect rejects missing columns", {
    bad_query <- data.frame(va = "TRAV1-1*01", cdr3a = "CAVRDSSYKLIF",
                            stringsAsFactors = FALSE)
    ref <- data.frame(
        va = "TRAV1-1*01", cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV19*01",  cdr3b = "CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    expect_error(tcrdist_rect(bad_query, ref, "human"), "missing required columns")
})

test_that("tcrdist_rect rejects NA values", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", NA_character_),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )
    ref <- tcrs[1L, ]
    expect_error(tcrdist_rect(tcrs, ref, "human"), "NA values")
})

test_that("tcrdist_rect returns correct row and column names", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    rect <- tcrdist_rect(tcrs[1L, ], tcrs, "human")
    expect_equal(rownames(rect), "1")
    expect_equal(colnames(rect), c("1", "2"))
})
