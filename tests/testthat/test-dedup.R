test_that("deduplicate=TRUE merges within-subject duplicates and sums counts", {
    tcrs <- data.frame(
        va      = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a   = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb      = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b   = c("CASSIRSSYEQYF", "CASSIRSSYEQYF", "CSARDRTGNTIYF"),
        subject = c("S1", "S1", "S1"),
        count   = c(3L, 5L, 2L),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", deduplicate = TRUE)
    df <- obj@clone_df
    expect_equal(nrow(df), 2L)

    # Merged row should have summed count
    merged <- df[df$cdr3b == "CASSIRSSYEQYF", ]
    expect_equal(nrow(merged), 1L)
    expect_equal(merged$count, 8L)

    # Non-duplicate row unchanged
    other <- df[df$cdr3b == "CSARDRTGNTIYF", ]
    expect_equal(other$count, 2L)
})

test_that("deduplicate=TRUE keeps cross-subject duplicates separate", {
    tcrs <- data.frame(
        va      = c("TRAV1-1*01", "TRAV1-1*01"),
        cdr3a   = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb      = c("TRBV19*01", "TRBV19*01"),
        cdr3b   = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        subject = c("S1", "S2"),
        count   = c(3L, 5L),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", deduplicate = TRUE)
    expect_equal(nrow(obj@clone_df), 2L)
})

test_that("deduplicate=TRUE without subject column merges all identical TCRs", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSSYEQYF", "CSARDRTGNTIYF"),
        count = c(3L, 5L, 2L),
        stringsAsFactors = FALSE
    )

    # No subject column -> groups by chain cols only
    obj <- TCRrep(tcrs, "human", deduplicate = TRUE)
    df <- obj@clone_df
    expect_equal(nrow(df), 2L)
    expect_equal(df$count[df$cdr3b == "CASSIRSSYEQYF"], 8L)
})

test_that("deduplicate=FALSE preserves all rows unchanged", {
    tcrs <- data.frame(
        va      = c("TRAV1-1*01", "TRAV1-1*01"),
        cdr3a   = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb      = c("TRBV19*01", "TRBV19*01"),
        cdr3b   = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        subject = c("S1", "S1"),
        count   = c(3L, 5L),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", deduplicate = FALSE)
    expect_equal(nrow(obj@clone_df), 2L)
    expect_equal(obj@clone_df$count, c(3L, 5L))
})

test_that("deduplicate with custom additional grouping columns", {
    tcrs <- data.frame(
        va      = c("TRAV1-1*01", "TRAV1-1*01", "TRAV1-1*01"),
        cdr3a   = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb      = c("TRBV19*01", "TRBV19*01", "TRBV19*01"),
        cdr3b   = c("CASSIRSSYEQYF", "CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        subject = c("S1", "S2", "S1"),
        count   = c(3L, 5L, 2L),
        stringsAsFactors = FALSE
    )

    # character(0): chain cols only, ignore subject -> all 3 merge into 1
    obj <- TCRrep(tcrs, "human", deduplicate = character(0))
    expect_equal(nrow(obj@clone_df), 1L)
    expect_equal(obj@clone_df$count, 10L)

    # c("subject"): chain cols + subject -> S1 rows merge, S2 stays separate
    obj2 <- TCRrep(tcrs, "human", deduplicate = c("subject"))
    expect_equal(nrow(obj2@clone_df), 2L)
})

test_that("deduplicate adds count=1 when count column is missing", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", deduplicate = TRUE)
    df <- obj@clone_df
    expect_equal(nrow(df), 1L)
    expect_true("count" %in% colnames(df))
    expect_equal(df$count, 2L)
})

test_that("deduplicate drops clone_id column after merge", {
    tcrs <- data.frame(
        va       = c("TRAV1-1*01", "TRAV1-1*01"),
        cdr3a    = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb       = c("TRBV19*01", "TRBV19*01"),
        cdr3b    = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        clone_id = c("clone_001", "clone_002"),
        count    = c(1L, 1L),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", deduplicate = TRUE)
    expect_false("clone_id" %in% colnames(obj@clone_df))
})

test_that("deduplicate drops NA rows with message", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", NA, "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    expect_message(
        obj <- TCRrep(tcrs, "human", deduplicate = TRUE),
        "dropping 1 row"
    )
    expect_equal(nrow(obj@clone_df), 2L)
})

test_that("deduplicate with no actual duplicates returns unchanged data", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CSARDRTGNTIYF"),
        count = c(3L, 5L),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", deduplicate = TRUE)
    expect_equal(nrow(obj@clone_df), 2L)
    expect_equal(obj@clone_df$count, c(3L, 5L))
})

test_that("deduplicate errors on nonexistent custom columns", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF"),
        vb    = c("TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )

    expect_error(
        TCRrep(tcrs, "human", deduplicate = c("nonexistent")),
        "not found"
    )
})

test_that("deduplicate works with alpha-only chains", {
    tcrs <- data.frame(
        va      = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a   = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        subject = c("S1", "S1", "S1"),
        count   = c(2L, 3L, 1L),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", chains = "A", deduplicate = TRUE)
    expect_equal(nrow(obj@clone_df), 2L)
    expect_equal(obj@clone_df$count[obj@clone_df$cdr3a == "CAVRDSSYKLIF"], 5L)
})

test_that("deduplicate on DASH dataset produces expected counts", {
    data(dash)

    # Default: chain cols + subject -> merges within-subject dups
    obj1 <- TCRrep(dash, "mouse", deduplicate = TRUE)
    expect_equal(nrow(obj1@clone_df), 1888L)

    # Chain cols only (no extra columns) -> merges cross-subject too
    obj2 <- TCRrep(dash, "mouse", deduplicate = character(0))
    expect_equal(nrow(obj2@clone_df), 1746L)

    # No dedup
    obj3 <- TCRrep(dash, "mouse", deduplicate = FALSE)
    expect_equal(nrow(obj3@clone_df), 1924L)
})
