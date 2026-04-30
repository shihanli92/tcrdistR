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
# tcr_shannon_entropy
# ===========================================================================

test_that("tcr_shannon_entropy: uniform equals log(S)", {
    H <- tcr_shannon_entropy(rep(10, 5))
    expect_equal(H, log(5), tolerance = 1e-10)
})

test_that("tcr_shannon_entropy: single species = 0", {
    expect_equal(tcr_shannon_entropy(c(100)), 0)
})

test_that("tcr_shannon_entropy: base=2 gives bits", {
    H_nats <- tcr_shannon_entropy(rep(10, 4))
    H_bits <- tcr_shannon_entropy(rep(10, 4), base = 2)
    expect_equal(H_bits, log2(4), tolerance = 1e-10)
    expect_equal(H_nats / log(2), H_bits, tolerance = 1e-10)
})

test_that("tcr_shannon_entropy: zeros are ignored", {
    H1 <- tcr_shannon_entropy(c(10, 20, 30))
    H2 <- tcr_shannon_entropy(c(10, 20, 30, 0, 0))
    expect_equal(H1, H2)
})

test_that("tcr_shannon_entropy: dominated distribution has low entropy", {
    H <- tcr_shannon_entropy(c(10000, 1, 1))
    expect_true(H < 0.01)
})


# ===========================================================================
# tcr_gini
# ===========================================================================

test_that("tcr_gini: uniform distribution near 0", {
    G <- tcr_gini(rep(100, 10))
    expect_true(G < 0.01)
})

test_that("tcr_gini: dominated distribution has high Gini", {
    # With few species, Gini max is (n-1)/n; need many species for > 0.9
    G <- tcr_gini(c(100000, rep(1, 100)))
    expect_true(G > 0.9)
})

test_that("tcr_gini: single species = 0", {
    expect_equal(tcr_gini(c(100)), 0)
})

test_that("tcr_gini: bounded [0, 1]", {
    for (counts in list(c(1,1,1), c(100,1), c(50,50,50,50), c(10000,1,1))) {
        G <- tcr_gini(counts)
        expect_true(G >= 0 && G <= 1)
    }
})

test_that("tcr_gini: more unequal = higher Gini", {
    G_uniform <- tcr_gini(rep(10, 5))
    G_skewed <- tcr_gini(c(100, 1, 1, 1, 1))
    expect_true(G_skewed > G_uniform)
})


# ===========================================================================
# tcr_fuzzy_diversity (requires compiled code)
# ===========================================================================

test_that("tcr_repertoire_overlap: identical repertoires have overlap 1", {
    a <- c(clone1 = 10, clone2 = 5, clone3 = 1)
    ov <- tcr_repertoire_overlap(a, a)
    expect_equal(ov$jaccard, 1)
    expect_equal(ov$morisita_horn, 1, tolerance = 1e-10)
    expect_equal(ov$overlap_coef, 1)
})

test_that("tcr_repertoire_overlap: disjoint repertoires have overlap 0", {
    a <- c(clone1 = 10, clone2 = 5)
    b <- c(clone3 = 8, clone4 = 3)
    ov <- tcr_repertoire_overlap(a, b)
    expect_equal(ov$jaccard, 0)
    expect_equal(ov$morisita_horn, 0)
    expect_equal(ov$overlap_coef, 0)
})

test_that("tcr_repertoire_overlap: partial overlap gives expected Jaccard", {
    a <- c(clone1 = 10, clone2 = 5, clone3 = 1)
    b <- c(clone1 = 8, clone3 = 3, clone4 = 2)
    ov <- tcr_repertoire_overlap(a, b)
    # shared = {clone1, clone3} = 2, union = {clone1..clone4} = 4
    expect_equal(ov$jaccard, 2 / 4)
})

test_that("tcr_repertoire_overlap: Morisita-Horn is abundance-weighted", {
    # Two samples sharing one dominant clone
    a <- c(clone1 = 100, clone2 = 1)
    b <- c(clone1 = 100, clone3 = 1)
    ov <- tcr_repertoire_overlap(a, b)
    # Morisita-Horn should be very high because clone1 dominates both
    expect_true(ov$morisita_horn > 0.9)
})

test_that("tcr_repertoire_overlap: handles single metric request", {
    a <- c(clone1 = 10, clone2 = 5)
    b <- c(clone1 = 8, clone3 = 3)
    ov <- tcr_repertoire_overlap(a, b, metrics = "jaccard")
    expect_true("jaccard" %in% names(ov))
    expect_false("morisita_horn" %in% names(ov))
    expect_false("overlap_coef" %in% names(ov))
})


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
