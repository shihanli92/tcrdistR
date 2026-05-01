test_that("tcrdist_sparse with Inf threshold matches tcrdist_matrix", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    sp  <- tcrdist_sparse(tcrs, "human", threshold = Inf)
    mat <- tcrdist_matrix(tcrs, "human")
    # Compare non-diagonal entries (diagonal is 0 in mat, structurally absent in sp)
    sp_dense <- as.matrix(sp)
    diag(mat)      <- 0
    diag(sp_dense) <- 0
    expect_equal(unname(sp_dense), unname(mat), tolerance = 1e-10)
})

test_that("tcrdist_sparse with threshold=100 captures only the close pair", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    # dist(1,2)=12, dist(1,3)=285, dist(2,3)=288 from tcrdist_matrix
    # Only pair (1,2) with dist=12 is within threshold=100
    sp <- tcrdist_sparse(tcrs, "human", threshold = 100)
    sp_dense <- as.matrix(sp)
    expect_equal(sp_dense[1, 2], 12)
    expect_equal(sp_dense[2, 1], 12)
    expect_equal(sp_dense[1, 3], 0)
    expect_equal(sp_dense[3, 1], 0)
    expect_equal(sp_dense[2, 3], 0)
    expect_equal(sp_dense[3, 2], 0)
})

test_that("tcrdist_sparse with threshold=250 captures all pairs", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    sp <- tcrdist_sparse(tcrs, "human", threshold = 300)
    sp_dense <- as.matrix(sp)
    # All three off-diagonal pairs should be present (all distances <= 300)
    mat <- tcrdist_matrix(tcrs, "human")
    expect_equal(sp_dense[1, 2], mat[1, 2])
    expect_equal(sp_dense[1, 3], mat[1, 3])
    expect_equal(sp_dense[2, 3], mat[2, 3])
})

test_that("tcrdist_sparse returns a dgCMatrix", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    sp <- tcrdist_sparse(tcrs, "human", threshold = 50)
    expect_true(inherits(sp, "dgCMatrix"))
})

test_that("tcrdist_sparse handles empty input", {
    tcrs <- data.frame(
        va    = character(0),
        cdr3a = character(0),
        vb    = character(0),
        cdr3b = character(0),
        stringsAsFactors = FALSE
    )
    sp <- tcrdist_sparse(tcrs, "human", threshold = 50)
    expect_true(inherits(sp, "dgCMatrix"))
    expect_equal(dim(sp), c(0L, 0L))
})

test_that("tcrdist_sparse rejects negative threshold", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF"),
        vb    = c("TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )
    expect_error(tcrdist_sparse(tcrs, "human", threshold = -1), "non-negative")
})

test_that("tcrdist_sparse threshold=0 returns only identical TCR pairs", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )
    # Identical TCRs have distance 0, so threshold=0 should capture them
    sp <- tcrdist_sparse(tcrs, "human", threshold = 0)
    sp_dense <- as.matrix(sp)
    expect_equal(sp_dense[1, 2], 0)
    expect_equal(sp_dense[2, 1], 0)
})
