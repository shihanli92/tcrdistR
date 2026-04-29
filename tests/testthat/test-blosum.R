test_that("BSD4 matrix has correct dimensions and names", {
    bsd4 <- bsd4_matrix()
    expect_equal(nrow(bsd4), 20L)
    expect_equal(ncol(bsd4), 20L)
    expect_equal(rownames(bsd4), AMINO_ACIDS)
    expect_equal(colnames(bsd4), AMINO_ACIDS)
})

test_that("BSD4 diagonal is all zeros", {
    bsd4 <- bsd4_matrix()
    expect_equal(unname(diag(bsd4)), rep(0.0, 20))
})

test_that("BSD4 is symmetric", {
    bsd4 <- bsd4_matrix()
    expect_true(isSymmetric(bsd4))
})

test_that("BSD4 values match BLOSUM62 derivation", {
    bsd4 <- bsd4_matrix()
    # A,A: BLOSUM62=4, diagonal -> 0
    expect_equal(bsd4["A", "A"], 0.0)
    # A,C: BLOSUM62=0, off-diag >=0 -> 4-0=4
    expect_equal(bsd4["A", "C"], 4.0)
    # F,Y: BLOSUM62=3, off-diag >=0 -> 4-3=1
    expect_equal(bsd4["F", "Y"], 1.0)
    # I,V: BLOSUM62=3, off-diag >=0 -> 4-3=1
    expect_equal(bsd4["I", "V"], 1.0)
    # W,W: BLOSUM62=11, diagonal -> 0
    expect_equal(bsd4["W", "W"], 0.0)
    # D,E: BLOSUM62=2, off-diag >=0 -> 4-2=2
    expect_equal(bsd4["D", "E"], 2.0)
    # A,W: BLOSUM62=-3, off-diag <0 -> 4
    expect_equal(bsd4["A", "W"], 4.0)
})

test_that("BSD4 values are non-negative and bounded", {
    bsd4 <- bsd4_matrix()
    expect_true(all(bsd4 >= 0))
    expect_true(all(bsd4 <= 4))
})

test_that("BLOSUM62 lookup returns correct values", {
    expect_equal(rcpp_blosum62_lookup("A", "A"), 4L)
    expect_equal(rcpp_blosum62_lookup("W", "W"), 11L)
    expect_equal(rcpp_blosum62_lookup("F", "Y"), 3L)
    expect_equal(rcpp_blosum62_lookup("A", "C"), 0L)
})

test_that("BLOSUM62 lookup is symmetric", {
    expect_equal(rcpp_blosum62_lookup("A", "W"), rcpp_blosum62_lookup("W", "A"))
    expect_equal(rcpp_blosum62_lookup("D", "E"), rcpp_blosum62_lookup("E", "D"))
})

test_that("BLOSUM62 lookup rejects invalid input", {
    expect_error(rcpp_blosum62_lookup("A", "Z"))
    expect_error(rcpp_blosum62_lookup("AB", "C"))
    expect_error(rcpp_blosum62_lookup("", "A"))
})
