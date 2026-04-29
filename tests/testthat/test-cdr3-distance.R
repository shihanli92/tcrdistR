test_that("identical CDR3 sequences have distance 0", {
    d <- weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSSYEQYF")
    expect_equal(d, 0.0)
})

test_that("CDR3 distance is symmetric", {
    d1 <- weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSYEQYF")
    d2 <- weighted_cdr3_distance("CASSIRSYEQYF", "CASSIRSSYEQYF")
    expect_equal(d1, d2)
})

test_that("CDR3 distance includes length penalty", {
    # Sequences differ by 1 character in length -> 12 gap penalty
    d <- weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSYEQYF")
    # d should include 1 * GAP_PENALTY_CDR3_REGION = 12
    expect_true(d >= 12)
})

test_that("CDR3 distance with custom weights", {
    # Same-length sequences so only alignment distance matters (no gap penalty)
    d_default <- weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRASYEQYF")
    d_w1 <- weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRASYEQYF",
                                    weight = 1L, gap_penalty = 12L)
    # With weight=1 instead of 3, alignment distance is 1/3 of default
    expect_true(d_w1 < d_default)
    expect_equal(d_w1 * 3L, d_default)
})

test_that("CDR3 distance rejects invalid input", {
    expect_error(weighted_cdr3_distance(NA_character_, "CASS"))
    expect_error(weighted_cdr3_distance("CASS", ""))
    expect_error(weighted_cdr3_distance(123, "CASS"))
    expect_error(weighted_cdr3_distance(c("CASS", "CASS"), "CASS"))
})

test_that("CDR3 distance with same-length sequences", {
    # Same length, different residues -> distance from BSD4 only, no gap penalty
    d <- weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRASYEQYF")
    # No length difference, so no gap penalty component
    # The substitution S->A at one position contributes
    expect_true(d > 0)
    expect_true(d < 12)  # less than one gap penalty
})

test_that("short CDR3 sequences get sentinel distance", {
    # Sequences with length < 5 after trimming should get sentinel
    d <- weighted_cdr3_distance("CAAA", "CAAAA")
    expect_true(d > 0)
})

test_that("rcpp_weighted_cdr3_distance matches R wrapper", {
    d_r <- weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSYEQYF")
    d_cpp <- rcpp_weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSYEQYF")
    expect_equal(d_r, d_cpp)
})
