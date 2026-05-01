# Tests for diversity metrics


# ===========================================================================
# tcr_diversity
# ===========================================================================

test_that("tcr_diversity computes correct entropy across distributions", {
    # Uniform: high entropy
    r_uniform <- tcr_diversity(rep(10, 5), order = 2)
    expect_true(r_uniform$entropy > 0.7 && r_uniform$entropy <= 1.0)
    expect_true(r_uniform$effective_number >= 1)

    # Single species: zero entropy
    r_single <- tcr_diversity(c(100), order = 2)
    expect_equal(r_single$entropy, 0)
    expect_equal(r_single$effective_number, 1)

    # Dominated: low entropy
    r_dom <- tcr_diversity(c(1000, 1, 1, 1), order = 2)
    expect_true(r_dom$entropy < 0.01)

    # Zeros are ignored
    r1 <- tcr_diversity(c(10, 20, 30), order = 2, ci = FALSE)
    r2 <- tcr_diversity(c(10, 20, 30, 0, 0), order = 2, ci = FALSE)
    expect_equal(r1$entropy, r2$entropy)

    # Different orders both yield valid results
    counts <- c(100, 10, 1, 1)
    r2o <- tcr_diversity(counts, order = 2, ci = FALSE)
    r3o <- tcr_diversity(counts, order = 3, ci = FALSE)
    expect_true(r2o$entropy >= 0 && r2o$entropy <= 1)
    expect_true(r3o$entropy >= 0 && r3o$entropy <= 1)

    # Confidence intervals bracket entropy
    r_ci <- tcr_diversity(c(50, 30, 20, 10, 5), order = 2, ci = TRUE)
    expect_true(r_ci$ci_lower <= r_ci$entropy)
    expect_true(r_ci$ci_upper >= r_ci$entropy)
    expect_true(r_ci$ci_lower >= 0 && r_ci$ci_upper <= 1)
})


# ===========================================================================
# tcr_richness, tcr_clonality, tcr_shannon_entropy, tcr_gini
# ===========================================================================

test_that("tcr_richness counts unique non-zero species", {
    expect_equal(tcr_richness(c(10, 5, 3, 1, 1)), 5L)
    expect_equal(tcr_richness(c(10, 0, 0, 5)), 2L)
    expect_equal(tcr_richness(c(100)), 1L)
})

test_that("tcr_clonality measures dominance correctly", {
    expect_true(tcr_clonality(rep(100, 10)) < 0.01)   # uniform ~ 0
    expect_true(tcr_clonality(c(10000, 1, 1)) > 0.9)  # dominated ~ 1
    expect_equal(tcr_clonality(c(100)), 1)              # single = 1
    # Always bounded [0, 1]
    for (counts in list(c(1, 1, 1), c(100, 1), c(50, 50, 50, 50))) {
        expect_true(tcr_clonality(counts) >= 0 && tcr_clonality(counts) <= 1)
    }
})

test_that("tcr_shannon_entropy has correct properties", {
    # Uniform = log(S)
    expect_equal(tcr_shannon_entropy(rep(10, 5)), log(5), tolerance = 1e-10)
    # Single species = 0
    expect_equal(tcr_shannon_entropy(c(100)), 0)
    # Base=2 gives bits
    H_bits <- tcr_shannon_entropy(rep(10, 4), base = 2)
    expect_equal(H_bits, log2(4), tolerance = 1e-10)
    # Zeros are ignored
    expect_equal(tcr_shannon_entropy(c(10, 20, 30)),
                 tcr_shannon_entropy(c(10, 20, 30, 0, 0)))
    # Dominated has low entropy
    expect_true(tcr_shannon_entropy(c(10000, 1, 1)) < 0.01)
})

test_that("tcr_gini measures inequality correctly", {
    expect_true(tcr_gini(rep(100, 10)) < 0.01)              # uniform ~ 0
    expect_true(tcr_gini(c(100000, rep(1, 100))) > 0.9)     # dominated ~ 1
    expect_equal(tcr_gini(c(100)), 0)                         # single = 0
    # Bounded [0, 1]
    for (counts in list(c(1, 1, 1), c(100, 1), c(50, 50, 50, 50))) {
        expect_true(tcr_gini(counts) >= 0 && tcr_gini(counts) <= 1)
    }
    # More unequal = higher Gini
    expect_true(tcr_gini(c(100, 1, 1, 1, 1)) > tcr_gini(rep(10, 5)))
})


# ===========================================================================
# tcr_repertoire_overlap
# ===========================================================================

test_that("tcr_repertoire_overlap handles identical, disjoint, and partial", {
    a <- c(clone1 = 10, clone2 = 5, clone3 = 1)
    # Identical
    ov_id <- tcr_repertoire_overlap(a, a)
    expect_equal(ov_id$jaccard, 1)
    expect_equal(ov_id$morisita_horn, 1, tolerance = 1e-10)
    expect_equal(ov_id$overlap_coef, 1)

    # Disjoint
    b_dis <- c(clone4 = 8, clone5 = 3)
    ov_dis <- tcr_repertoire_overlap(a, b_dis)
    expect_equal(ov_dis$jaccard, 0)
    expect_equal(ov_dis$morisita_horn, 0)

    # Partial overlap
    b_part <- c(clone1 = 8, clone3 = 3, clone4 = 2)
    ov_part <- tcr_repertoire_overlap(a, b_part)
    expect_equal(ov_part$jaccard, 2 / 4)

    # Morisita-Horn is abundance-weighted
    a_dom <- c(clone1 = 100, clone2 = 1)
    b_dom <- c(clone1 = 100, clone3 = 1)
    ov_dom <- tcr_repertoire_overlap(a_dom, b_dom)
    expect_true(ov_dom$morisita_horn > 0.9)

    # Single metric request
    ov_single <- tcr_repertoire_overlap(a, b_part, metrics = "jaccard")
    expect_true("jaccard" %in% names(ov_single))
    expect_false("morisita_horn" %in% names(ov_single))
})


# ===========================================================================
# tcr_fuzzy_diversity
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
    expect_true(result$fuzzy_diversity >= 0 && result$fuzzy_diversity <= 1)
    expect_true(result$fuzzy_diversity <= result$standard_diversity + 0.01)
})
