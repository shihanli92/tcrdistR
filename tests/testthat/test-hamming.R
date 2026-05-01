# Tests for Hamming distance functions


test_that("hamming_distance: identical sequences = 0", {
    expect_equal(hamming_distance("CASSI", "CASSI"), 0L)
})

test_that("hamming_distance: one mismatch = 1", {
    expect_equal(hamming_distance("CASSI", "CASSK"), 1L)
})

test_that("hamming_distance: all mismatches", {
    expect_equal(hamming_distance("AAA", "BBB"), 3L)
})

test_that("hamming_distance: different lengths = -1", {
    expect_equal(hamming_distance("CASSI", "CASSILY"), -1L)
})

test_that("hamming_distance: empty strings = 0", {
    expect_equal(hamming_distance("", ""), 0L)
})


test_that("hamming_matrix: symmetric", {
    seqs <- c("CASSI", "CASSK", "CASRL")
    mat <- hamming_matrix(seqs)
    expect_true(isSymmetric(mat))
})

test_that("hamming_matrix: diagonal is zero", {
    seqs <- c("CASSI", "CASSK", "CASRL")
    mat <- hamming_matrix(seqs)
    expect_true(all(diag(mat) == 0L))
})

test_that("hamming_matrix: consistent with single-pair", {
    seqs <- c("CASSI", "CASSK", "CASRL")
    mat <- hamming_matrix(seqs)
    for (i in seq_along(seqs)) {
        for (j in seq_along(seqs)) {
            d <- hamming_distance(seqs[i], seqs[j])
            if (d >= 0) {
                expect_equal(mat[i, j], d)
            }
        }
    }
})

test_that("hamming_matrix: different lengths get max penalty", {
    seqs <- c("CASSI", "CASSILY")  # 5 vs 7
    mat <- hamming_matrix(seqs)
    expect_equal(mat[1, 2], 7L)  # max(5, 7) = 7
    expect_equal(mat[2, 1], 7L)
})

test_that("hamming_matrix: single sequence", {
    mat <- hamming_matrix("CASSI")
    expect_equal(dim(mat), c(1L, 1L))
    expect_equal(mat[1, 1], 0L)
})
