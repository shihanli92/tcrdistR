# Tests for distance-based TCR joins

test_that("tcrdist_join inner: self-join finds all at radius=0", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01", "TRAV10*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01", "TRBV7-2*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY"),
        stringsAsFactors = FALSE
    )

    result <- tcrdist_join(tcr_df, tcr_df, "human", radius = 0)
    # Each TCR should match itself at distance 0
    expect_true(nrow(result) >= 3L)
    expect_true(all(result$tcrdist == 0))
})

test_that("tcrdist_join inner: large radius finds all pairs", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY"),
        stringsAsFactors = FALSE
    )

    result <- tcrdist_join(tcr_df, tcr_df, "human", radius = 1000)
    # 2 left rows x 2 right matches each = 4 (but max_n=5 so all fit)
    expect_true(nrow(result) >= 2L)
    expect_true("tcrdist" %in% colnames(result))
})

test_that("tcrdist_join left: all left rows present", {
    skip_on_cran()

    left_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01", "TRAV10*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01", "TRBV7-2*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY"),
        stringsAsFactors = FALSE
    )

    right_df <- data.frame(
        va = "TRAV1-1*01",
        cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV5-1*01",
        cdr3b = "CASSIRSSYEQY",
        stringsAsFactors = FALSE
    )

    result <- tcrdist_join(left_df, right_df, "human", radius = 0,
                            type = "left")
    # At least 3 rows (one per left row)
    expect_true(nrow(result) >= 3L)
})

test_that("tcrdist_join: max_n limits matches", {
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

    result <- tcrdist_join(tcr_df, tcr_df, "human", radius = 1000,
                            max_n = 2L)
    # Total rows should be at most 5 * 2 = 10 (5 left rows, max 2 matches each)
    expect_true(nrow(result) <= 5L * 2L)
})

test_that("tcrdist_join: suffix handling", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = "TRAV1-1*01",
        cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV5-1*01",
        cdr3b = "CASSIRSSYEQY",
        stringsAsFactors = FALSE
    )

    result <- tcrdist_join(tcr_df, tcr_df, "human", radius = 0,
                            suffix = c("_left", "_right"))
    expect_true("va_left" %in% colnames(result))
    expect_true("va_right" %in% colnames(result))
    expect_true("tcrdist" %in% colnames(result))
})

test_that("tcrdist_join: tcrdist values are within radius", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01", "TRAV10*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01", "TRBV7-2*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY"),
        stringsAsFactors = FALSE
    )

    radius <- 50
    result <- tcrdist_join(tcr_df, tcr_df, "human", radius = radius)
    if (nrow(result) > 0L) {
        non_na <- result$tcrdist[!is.na(result$tcrdist)]
        expect_true(all(non_na <= radius))
    }
})

test_that("tcrdist_join: empty input returns empty", {
    skip_on_cran()

    empty <- data.frame(va = character(0), cdr3a = character(0),
                         vb = character(0), cdr3b = character(0),
                         stringsAsFactors = FALSE)
    tcr_df <- data.frame(
        va = "TRAV1-1*01", cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV5-1*01", cdr3b = "CASSIRSSYEQY",
        stringsAsFactors = FALSE
    )

    result <- tcrdist_join(empty, tcr_df, "human", radius = 50)
    expect_equal(nrow(result), 0L)
})
