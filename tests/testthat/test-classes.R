test_that("TCRrep basic creation succeeds and returns valid object", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human")
    expect_true(is(obj, "TCRrep"))
    expect_equal(nrow(obj@clone_df), 3L)
    expect_equal(obj@organism, "human")
    expect_equal(obj@chains,   "AB")
    expect_equal(obj@metric,   "tcrdist")
})

test_that("TCRrep rejects invalid chains argument", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF"),
        vb    = c("TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )

    expect_error(TCRrep(tcrs, "human", chains = "XYZ"))
})

test_that("TCRrep rejects clone_df with missing required columns for chains='AB'", {
    # Only alpha columns — missing vb and cdr3b
    tcrs_alpha_only <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        stringsAsFactors = FALSE
    )

    expect_error(TCRrep(tcrs_alpha_only, "human", chains = "AB"))
})

test_that("TCRrep single-chain 'A' works with only alpha columns", {
    tcrs_alpha <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs_alpha, "human", chains = "A", deduplicate = FALSE)
    expect_true(is(obj, "TCRrep"))
    expect_equal(obj@chains, "A")
    expect_equal(nrow(obj@clone_df), 3L)
})

test_that("TCRrep coerces factor columns to character", {
    tcrs_factor <- data.frame(
        va    = factor(c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01")),
        cdr3a = factor(c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF")),
        vb    = factor(c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01")),
        cdr3b = factor(c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"))
    )

    obj <- TCRrep(tcrs_factor, "human")
    expect_true(is(obj, "TCRrep"))
    expect_true(is.character(obj@clone_df$va))
    expect_true(is.character(obj@clone_df$cdr3a))
    expect_true(is.character(obj@clone_df$vb))
    expect_true(is.character(obj@clone_df$cdr3b))
})

test_that("TCRrep with compute_distances=TRUE stores a matrix in paired_dist", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", compute_distances = TRUE)
    expect_false(is.null(obj@paired_dist))
    expect_true(is.matrix(obj@paired_dist))
    expect_equal(dim(obj@paired_dist), c(3L, 3L))
})

test_that("TCRrep with compute_distances=FALSE (default) leaves paired_dist NULL", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human")
    expect_true(is.null(obj@paired_dist))
})

test_that("show method produces expected output for TCRrep", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human")
    out <- capture.output(show(obj))
    expect_true(any(grepl("TCRrep object", out)))
    expect_true(any(grepl("clonotypes", out)))
})

test_that("TCRrep succeeds with empty data.frame", {
    empty_df <- data.frame(
        va    = character(0L),
        cdr3a = character(0L),
        vb    = character(0L),
        cdr3b = character(0L),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(empty_df, "human")
    expect_true(is(obj, "TCRrep"))
    expect_equal(nrow(obj@clone_df), 0L)
})

test_that("TCRrep default weights and gap penalties match package constants", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF"),
        vb    = c("TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human")
    expect_equal(obj@weights$cdr3,          3L)
    expect_equal(obj@weights$v_region,       1L)
    expect_equal(obj@gap_penalties$cdr3,    12L)
    expect_equal(obj@gap_penalties$v_region,  4L)
})

test_that("TCRrep rejects removed metrics (levenshtein, nw)", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF"),
        vb    = c("TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )

    expect_error(TCRrep(tcrs, "human", metric = "levenshtein"))
    expect_error(TCRrep(tcrs, "human", metric = "nw"))
})

test_that("TCRrep with metric='hamming' computes hamming distances", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    obj <- TCRrep(tcrs, "human", metric = "hamming", compute_distances = TRUE)
    expect_true(is(obj, "TCRrep"))
    expect_equal(obj@metric, "hamming")
    expect_false(is.null(obj@paired_dist))
    expect_true(is.matrix(obj@paired_dist))
    expect_equal(dim(obj@paired_dist), c(3L, 3L))
    # Diagonal should be 0
    expect_equal(diag(obj@paired_dist), c(0, 0, 0))
})


# ===========================================================================
# Subsetting: [ and subset()
# ===========================================================================

test_that("[ subsets clone_df and paired_dist by integer indices", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:20, ], "mouse", compute_distances = TRUE)
    sub <- rep[1:5, ]
    expect_equal(nrow(sub@clone_df), 5L)
    expect_equal(dim(sub@paired_dist), c(5L, 5L))
    expect_equal(sub@paired_dist, rep@paired_dist[1:5, 1:5])
})

test_that("[ subsets by logical vector", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:20, ], "mouse", compute_distances = TRUE)
    idx <- rep@clone_df$epitope == rep@clone_df$epitope[1]
    sub <- rep[idx, ]
    expect_equal(nrow(sub@clone_df), sum(idx))
    expect_true(all(sub@clone_df$epitope == rep@clone_df$epitope[1]))
    expect_equal(dim(sub@paired_dist), c(sum(idx), sum(idx)))
})

test_that("[ rejects logical vector of wrong length", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:10, ], "mouse")
    expect_error(rep[c(TRUE, FALSE), ], "logical index length")
})

test_that("[ clears knn and meta_clonotypes with message", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:10, ], "mouse", compute_distances = TRUE)
    rep@knn_indices <- matrix(1L, nrow = 10, ncol = 2)
    rep@knn_distances <- matrix(0, nrow = 10, ncol = 2)
    rep@meta_clonotypes <- data.frame(x = 1)
    expect_message(
        expect_message(sub <- rep[1:5, ], "KNN"),
        "meta-clonotypes"
    )
    expect_null(sub@knn_indices)
    expect_null(sub@knn_distances)
    expect_null(sub@meta_clonotypes)
})

test_that("[ does not message when knn/meta are already NULL", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:10, ], "mouse", compute_distances = TRUE)
    expect_no_message(sub <- rep[1:5, ])
})

test_that("[ preserves organism, chains, metric, weights", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:10, ], "mouse", compute_distances = TRUE)
    sub <- rep[1:3, ]
    expect_equal(sub@organism, "mouse")
    expect_equal(sub@chains, "AB")
    expect_equal(sub@metric, "tcrdist")
    expect_equal(sub@weights, rep@weights)
    expect_equal(sub@gap_penalties, rep@gap_penalties)
})

test_that("[ works without distances computed", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:10, ], "mouse", compute_distances = FALSE)
    sub <- rep[1:5, ]
    expect_equal(nrow(sub@clone_df), 5L)
    expect_null(sub@paired_dist)
})

test_that("subset() filters by expression", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:50, ], "mouse", compute_distances = TRUE)
    epi <- rep@clone_df$epitope[1]
    sub <- subset(rep, epitope == epi)
    expect_true(all(sub@clone_df$epitope == epi))
    n <- sum(rep@clone_df$epitope == epi)
    expect_equal(nrow(sub@clone_df), n)
    expect_equal(dim(sub@paired_dist), c(n, n))
})

test_that("subset() distances match [ subsetting", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:50, ], "mouse", compute_distances = TRUE)
    epi <- rep@clone_df$epitope[1]
    sub1 <- subset(rep, epitope == epi)
    idx <- which(rep@clone_df$epitope == epi)
    sub2 <- rep[idx, ]
    expect_equal(sub1@paired_dist, sub2@paired_dist)
})

test_that("subset() handles NA in filter column", {
    data(dash, envir = environment())
    tcr_rep <- TCRrep(dash[1:10, ], "mouse")
    tcr_rep@clone_df$epitope[1] <- NA
    target_epi <- tcr_rep@clone_df$epitope[2]
    sub <- subset(tcr_rep, epitope == target_epi)
    expect_false(any(is.na(sub@clone_df$epitope)))
})


# ===========================================================================
# Gamma-delta (GD) chain mode
# ===========================================================================

test_that("TCRrep with chains='GD' and organism='human_gd' creates valid object", {
    gd_df <- data.frame(
        va    = c("TRGV1*01",  "TRGV1*01",  "TRGV10*01"),
        cdr3a = c("CATWDRF",   "CATWDRF",   "CATWDSF"),
        vb    = c("TRDV1*01",  "TRDV2*01",  "TRDV1*01"),
        cdr3b = c("CALGELGDDKLIF", "CALGELGDDKLIF", "CALGELSDDKLIF"),
        stringsAsFactors = FALSE
    )
    obj <- TCRrep(gd_df, "human_gd", chains = "GD")
    expect_true(is(obj, "TCRrep"))
    expect_equal(obj@chains, "GD")
    expect_equal(obj@organism, "human_gd")
})

test_that("TCRrep GD with compute_distances=TRUE produces distance matrix", {
    gd_df <- data.frame(
        va    = c("TRGV1*01",  "TRGV1*01",  "TRGV10*01"),
        cdr3a = c("CATWDRF",   "CATWDRF",   "CATWDSF"),
        vb    = c("TRDV1*01",  "TRDV2*01",  "TRDV1*01"),
        cdr3b = c("CALGELGDDKLIF", "CALGELGDDKLIF", "CALGELSDDKLIF"),
        stringsAsFactors = FALSE
    )
    obj <- TCRrep(gd_df, "human_gd", chains = "GD", compute_distances = TRUE)
    expect_false(is.null(obj@paired_dist))
    expect_true(is.matrix(obj@paired_dist))
    expect_equal(dim(obj@paired_dist), c(3L, 3L))
    expect_equal(unname(diag(obj@paired_dist)), c(0, 0, 0))
})

test_that("TCRrep GD with mouse_gd organism works", {
    gd_df <- data.frame(
        va    = c("TRGV1*01",  "TRGV1*01"),
        cdr3a = c("CATWDRF",   "CATWDSF"),
        vb    = c("TRAV1*01",  "TRAV1*01"),
        cdr3b = c("CALGELGDDKLIF", "CALGELSDDKLIF"),
        stringsAsFactors = FALSE
    )
    obj <- TCRrep(gd_df, "mouse_gd", chains = "GD")
    expect_true(is(obj, "TCRrep"))
    expect_equal(obj@organism, "mouse_gd")
})
