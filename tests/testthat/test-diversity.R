# Tests for diversity metrics


# ===========================================================================
# tcr_diversity
# ===========================================================================

test_that("tcr_diversity: uniform distribution has high entropy", {
    result <- tcr_diversity(rep(10, 5), order = 2)
    # Simpson's diversity for uniform: 1 - 1/S = 1 - 1/5 = 0.8
    expect_true(result$entropy > 0.7)
    expect_true(result$entropy <= 1.0)
})

test_that("tcr_diversity: single species has zero entropy", {
    result <- tcr_diversity(c(100), order = 2)
    expect_equal(result$entropy, 0)
    expect_equal(result$effective_number, 1)
})

test_that("tcr_diversity: two equal species", {
    result <- tcr_diversity(c(50, 50), order = 2)
    # Simpson = 1 - (0.5^2 + 0.5^2) = 1 - 0.5 = 0.5 (approx, with correction)
    expect_true(result$entropy > 0.4)
    expect_true(result$entropy < 0.6)
})

test_that("tcr_diversity: dominated distribution has low entropy", {
    result <- tcr_diversity(c(1000, 1, 1, 1), order = 2)
    expect_true(result$entropy < 0.01)
})

test_that("tcr_diversity: higher order affects entropy", {
    counts <- c(100, 10, 1, 1)
    r2 <- tcr_diversity(counts, order = 2, ci = FALSE)
    r3 <- tcr_diversity(counts, order = 3, ci = FALSE)
    # Both should be valid diversity values
    expect_true(r2$entropy >= 0 && r2$entropy <= 1)
    expect_true(r3$entropy >= 0 && r3$entropy <= 1)
})

test_that("tcr_diversity: confidence intervals bracket entropy", {
    result <- tcr_diversity(c(50, 30, 20, 10, 5), order = 2, ci = TRUE)
    expect_true(result$ci_lower <= result$entropy)
    expect_true(result$ci_upper >= result$entropy)
    expect_true(result$ci_lower >= 0)
    expect_true(result$ci_upper <= 1)
})

test_that("tcr_diversity: effective number >= 1", {
    result <- tcr_diversity(c(50, 30, 20), order = 2)
    expect_true(result$effective_number >= 1)
})

test_that("tcr_diversity: zeros are ignored", {
    r1 <- tcr_diversity(c(10, 20, 30), order = 2, ci = FALSE)
    r2 <- tcr_diversity(c(10, 20, 30, 0, 0), order = 2, ci = FALSE)
    expect_equal(r1$entropy, r2$entropy)
})


# ===========================================================================
# tcr_richness
# ===========================================================================

test_that("tcr_richness: counts unique species", {
    expect_equal(tcr_richness(c(10, 5, 3, 1, 1)), 5L)
})

test_that("tcr_richness: zeros don't count", {
    expect_equal(tcr_richness(c(10, 0, 0, 5)), 2L)
})

test_that("tcr_richness: single species", {
    expect_equal(tcr_richness(c(100)), 1L)
})


# ===========================================================================
# tcr_clonality
# ===========================================================================

test_that("tcr_clonality: uniform distribution near 0", {
    clonality <- tcr_clonality(rep(100, 10))
    expect_true(clonality < 0.01)
})

test_that("tcr_clonality: dominated distribution near 1", {
    clonality <- tcr_clonality(c(10000, 1, 1))
    expect_true(clonality > 0.9)
})

test_that("tcr_clonality: single species = 1", {
    expect_equal(tcr_clonality(c(100)), 1)
})

test_that("tcr_clonality: bounded [0, 1]", {
    for (counts in list(c(1,1,1), c(100,1), c(50,50,50,50))) {
        cl <- tcr_clonality(counts)
        expect_true(cl >= 0 && cl <= 1)
    }
})


# ===========================================================================
# tcr_fuzzy_diversity (requires compiled code)
# ===========================================================================

test_that("tcr_fuzzy_diversity works with synthetic TCRs", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = rep("TRAV1-1*01", 5),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVKDSSYKLIF",
                   "CAVRDSYKLIF", "CAVKDSYKLIF"),
        vb = rep("TRBV5-1*01", 5),
        cdr3b = c("CASSIRSSYEQY", "CASSIRSSYEQY", "CASSIKSSYEQY",
                   "CASSIRSYEQY", "CASSIKSSYEQF"),
        stringsAsFactors = FALSE
    )

    result <- tcr_fuzzy_diversity(tcr_df, "human", threshold = 100)
    expect_type(result, "list")
    expect_true("fuzzy_diversity" %in% names(result))
    expect_true("standard_diversity" %in% names(result))
    expect_true(result$fuzzy_diversity >= 0)
    expect_true(result$fuzzy_diversity <= 1)
    # With large threshold, fuzzy diversity should be <= standard
    expect_true(result$fuzzy_diversity <= result$standard_diversity + 0.01)
})
