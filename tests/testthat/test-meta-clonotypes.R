# Tests for meta-clonotype detection


# ===========================================================================
# .compute_distance_ecdf
# ===========================================================================

test_that(".compute_distance_ecdf: monotonically increasing", {
    ecdf_fn <- tcrdistR:::.compute_distance_ecdf(c(5, 10, 15, 20, 25))
    vals <- ecdf_fn(seq(0, 30, by = 5))
    expect_true(all(diff(vals) >= 0))
})

test_that(".compute_distance_ecdf: bounded [0, 1]", {
    ecdf_fn <- tcrdistR:::.compute_distance_ecdf(c(1, 2, 3, 4, 5))
    vals <- ecdf_fn(seq(0, 10, by = 1))
    expect_true(all(vals >= 0))
    expect_true(all(vals <= 1))
})

test_that(".compute_distance_ecdf: correct at extremes", {
    ecdf_fn <- tcrdistR:::.compute_distance_ecdf(c(5, 10, 15))
    expect_equal(ecdf_fn(0), 0)  # below all values
    expect_equal(ecdf_fn(100), 1)  # above all values
})

test_that(".compute_distance_ecdf: weighted ECDF", {
    ecdf_fn <- tcrdistR:::.compute_distance_ecdf(
        c(5, 10, 15), weights = c(1, 2, 1))
    # At d=10: (1+2)/4 = 0.75
    expect_equal(ecdf_fn(10), 0.75)
})

test_that(".compute_distance_ecdf: empty input", {
    ecdf_fn <- tcrdistR:::.compute_distance_ecdf(numeric(0))
    expect_equal(ecdf_fn(5), 0)
})


# ===========================================================================
# .calc_radii
# ===========================================================================

test_that(".calc_radii: returns integer vector of correct length", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY"),
        stringsAsFactors = FALSE
    )

    bg_df <- data.frame(
        va = rep("TRAV10*01", 20),
        cdr3a = paste0("CAVRD", sample(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"), 20, replace = TRUE), "YKLIF"),
        vb = rep("TRBV7-2*01", 20),
        cdr3b = paste0("CASSIR", sample(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"), 20, replace = TRUE), "YEQY"),
        stringsAsFactors = FALSE
    )

    radii <- tcrdistR:::.calc_radii(tcr_df, "human", bg_df,
                                     ctrl_bkgd = 0.5, max_radius = 200L)
    expect_length(radii, 2L)
    expect_true(is.integer(radii) || is.numeric(radii))
    expect_true(all(radii <= 200L))
    expect_true(all(radii >= 0L))
})


# ===========================================================================
# find_meta_clonotypes
# ===========================================================================

test_that("find_meta_clonotypes: returns correct structure", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = rep("TRAV1-1*01", 6),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVKDSSYKLIF",
                   "CAVRDSYKLIF", "CAVKDSYKLIF", "CAVRDSSYKLIF"),
        vb = rep("TRBV5-1*01", 6),
        cdr3b = c("CASSIRSSYEQY", "CASSIRSSYEQY", "CASSIKSSYEQY",
                   "CASSIRSYEQY", "CASSIKSSYEQF", "CASSIRSSYEQY"),
        subject = c("S1", "S2", "S1", "S2", "S1", "S3"),
        stringsAsFactors = FALSE
    )

    result <- find_meta_clonotypes(tcr_df, "human", radius = 100,
                                    min_nsubject = 2L)
    expect_s3_class(result, "data.frame")
    expected_cols <- c("center_index", "va", "cdr3a", "vb", "cdr3b",
                       "radius", "K_neighbors", "nsubject",
                       "neighbor_indices")
    expect_true(all(expected_cols %in% colnames(result)))
})

test_that("find_meta_clonotypes: filters by min_nsubject", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = rep("TRAV1-1*01", 4),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVKDSSYKLIF",
                   "CAVRDSYKLIF"),
        vb = rep("TRBV5-1*01", 4),
        cdr3b = c("CASSIRSSYEQY", "CASSIRSSYEQY", "CASSIKSSYEQY",
                   "CASSIRSYEQY"),
        subject = c("S1", "S1", "S1", "S1"),  # all same subject
        stringsAsFactors = FALSE
    )

    result <- find_meta_clonotypes(tcr_df, "human", radius = 100,
                                    min_nsubject = 2L)
    # No meta-clonotypes since all from same subject
    expect_equal(nrow(result), 0L)
})

test_that("find_meta_clonotypes: finds shared motifs across subjects", {
    skip_on_cran()

    # Two groups of identical TCRs from different subjects
    tcr_df <- data.frame(
        va = rep("TRAV1-1*01", 4),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF",
                   "CAVKDSSYKLIF", "CAVKDSSYKLIF"),
        vb = rep("TRBV5-1*01", 4),
        cdr3b = c("CASSIRSSYEQY", "CASSIRSSYEQY",
                   "CASSIKSSYEQY", "CASSIKSSYEQY"),
        subject = c("S1", "S2", "S1", "S2"),
        stringsAsFactors = FALSE
    )

    result <- find_meta_clonotypes(tcr_df, "human", radius = 0,
                                    min_nsubject = 2L)
    # Should find at least 1 meta-clonotype (identical TCRs from 2 subjects)
    expect_true(nrow(result) >= 1L)
    expect_true(all(result$nsubject >= 2L))
})

test_that("find_meta_clonotypes: errors on missing subject column", {
    tcr_df <- data.frame(
        va = "TRAV1-1*01", cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV5-1*01", cdr3b = "CASSIRSSYEQY",
        stringsAsFactors = FALSE
    )

    expect_error(
        find_meta_clonotypes(tcr_df, "human", radius = 50),
        "not found"
    )
})


# ===========================================================================
# summarize_meta_clonotype
# ===========================================================================

test_that("summarize_meta_clonotype returns correct structure", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-1*01", "TRAV1-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVKDSSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV5-1*01", "TRBV6-1*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIRSSYEQF", "CASSIKSSYEQY"),
        stringsAsFactors = FALSE
    )

    summary <- summarize_meta_clonotype(tcr_df, 1L, c(1L, 2L, 3L))
    expect_type(summary, "list")
    expect_equal(summary$n_members, 3L)
    expect_true("va_usage" %in% names(summary))
    expect_true("vb_usage" %in% names(summary))
    expect_length(summary$cdr3a_lengths, 3L)
})


# ===========================================================================
# dist_matrix bypass
# ===========================================================================

test_that("find_meta_clonotypes accepts precomputed dist_matrix", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-1*01", "TRAV1-2*01",
               "TRAV1-1*01", "TRAV1-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVKDSSYKLIF",
                   "CAVRDSSYKLIF", "CAVKDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV5-1*01", "TRBV6-1*01",
               "TRBV5-1*01", "TRBV6-1*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIRSSYEQF", "CASSIKSSYEQY",
                   "CASSIRSSYEQY", "CASSIKSSYEQF"),
        subject = c("S1", "S2", "S1", "S2", "S1"),
        stringsAsFactors = FALSE
    )

    dm <- tcrdist_matrix(tcr_df, "human")
    r1 <- find_meta_clonotypes(tcr_df, "human", radius = 200)
    r2 <- find_meta_clonotypes(tcr_df, "human", radius = 200,
                                dist_matrix = dm)
    expect_equal(r1, r2)
})
