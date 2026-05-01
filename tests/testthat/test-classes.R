# Tests for TCRrep S4 class


# ===========================================================================
# TCRrep creation, validation, and options
# ===========================================================================

test_that("TCRrep creation with various options", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01",  "TRAV1-1*01",  "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV19*01",   "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    # Basic creation
    obj <- TCRrep(tcrs, "human")
    expect_true(is(obj, "TCRrep"))
    expect_equal(nrow(obj@clone_df), 3L)
    expect_equal(obj@organism, "human")
    expect_equal(obj@chains, "AB")
    expect_equal(obj@metric, "tcrdist")
    expect_true(is.null(obj@paired_dist))

    # compute_distances=TRUE
    obj_dist <- TCRrep(tcrs, "human", compute_distances = TRUE)
    expect_false(is.null(obj_dist@paired_dist))
    expect_true(is.matrix(obj_dist@paired_dist))
    expect_equal(dim(obj_dist@paired_dist), c(3L, 3L))

    # Default weights and gap penalties
    expect_equal(obj@weights$cdr3, 3L)
    expect_equal(obj@weights$v_region, 1L)
    expect_equal(obj@gap_penalties$cdr3, 12L)
    expect_equal(obj@gap_penalties$v_region, 4L)

    # Factor coercion
    tcrs_fac <- tcrs
    tcrs_fac[] <- lapply(tcrs_fac, factor)
    obj_fac <- TCRrep(tcrs_fac, "human")
    expect_true(is.character(obj_fac@clone_df$va))
    expect_true(is.character(obj_fac@clone_df$cdr3b))

    # Empty data.frame
    empty_df <- data.frame(va = character(0), cdr3a = character(0),
                            vb = character(0), cdr3b = character(0),
                            stringsAsFactors = FALSE)
    obj_empty <- TCRrep(empty_df, "human")
    expect_equal(nrow(obj_empty@clone_df), 0L)

    # Single-chain A
    alpha_df <- tcrs[, c("va", "cdr3a")]
    obj_a <- TCRrep(alpha_df, "human", chains = "A", deduplicate = FALSE)
    expect_equal(obj_a@chains, "A")
    expect_equal(nrow(obj_a@clone_df), 3L)
})

test_that("TCRrep rejects invalid inputs", {
    tcrs <- data.frame(
        va = "TRAV1-1*01", cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV19*01", cdr3b = "CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    expect_error(TCRrep(tcrs, "human", chains = "XYZ"))
    expect_error(TCRrep(tcrs[, c("va", "cdr3a")], "human", chains = "AB"))
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
    expect_equal(obj@metric, "hamming")
    expect_false(is.null(obj@paired_dist))
    expect_equal(dim(obj@paired_dist), c(3L, 3L))
    expect_equal(diag(obj@paired_dist), c(0, 0, 0))
})


# ===========================================================================
# Subsetting: [ and subset()
# ===========================================================================

test_that("[ subsets correctly by integer and logical indices", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:20, ], "mouse", compute_distances = TRUE)

    # Integer subsetting
    sub_int <- rep[1:5, ]
    expect_equal(nrow(sub_int@clone_df), 5L)
    expect_equal(dim(sub_int@paired_dist), c(5L, 5L))
    expect_equal(sub_int@paired_dist, rep@paired_dist[1:5, 1:5])

    # Logical subsetting
    idx <- rep@clone_df$epitope == rep@clone_df$epitope[1]
    sub_log <- rep[idx, ]
    expect_equal(nrow(sub_log@clone_df), sum(idx))
    expect_true(all(sub_log@clone_df$epitope == rep@clone_df$epitope[1]))

    # Wrong-length logical rejected
    expect_error(rep[c(TRUE, FALSE), ], "logical index length")

    # Without distances
    rep_no_dist <- TCRrep(dash[1:10, ], "mouse", compute_distances = FALSE)
    sub_no_dist <- rep_no_dist[1:5, ]
    expect_equal(nrow(sub_no_dist@clone_df), 5L)
    expect_null(sub_no_dist@paired_dist)

    # Preserves organism, chains, metric, weights
    expect_equal(sub_int@organism, "mouse")
    expect_equal(sub_int@chains, "AB")
    expect_equal(sub_int@metric, "tcrdist")
    expect_equal(sub_int@weights, rep@weights)
    expect_equal(sub_int@gap_penalties, rep@gap_penalties)
})

test_that("[ clears knn/meta_clonotypes and messages appropriately", {
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
    expect_null(sub@meta_clonotypes)

    # No message when knn/meta already NULL
    rep2 <- TCRrep(dash[1:10, ], "mouse", compute_distances = TRUE)
    expect_no_message(rep2[1:5, ])
})

test_that("subset() filters and matches [ behavior", {
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:50, ], "mouse", compute_distances = TRUE)
    epi <- rep@clone_df$epitope[1]

    sub <- subset(rep, epitope == epi)
    expect_true(all(sub@clone_df$epitope == epi))
    n <- sum(rep@clone_df$epitope == epi)
    expect_equal(nrow(sub@clone_df), n)
    expect_equal(dim(sub@paired_dist), c(n, n))

    # Matches [ subsetting
    idx <- which(rep@clone_df$epitope == epi)
    sub2 <- rep[idx, ]
    expect_equal(sub@paired_dist, sub2@paired_dist)

    # NA handling
    rep_na <- TCRrep(dash[1:10, ], "mouse")
    rep_na@clone_df$epitope[1] <- NA
    target <- rep_na@clone_df$epitope[2]
    sub_na <- subset(rep_na, epitope == target)
    expect_false(any(is.na(sub_na@clone_df$epitope)))
})


# ===========================================================================
# Gamma-delta (GD) chain mode
# ===========================================================================

test_that("TCRrep GD mode works for human and mouse", {
    gd_df <- data.frame(
        va    = c("TRGV1*01",  "TRGV1*01",  "TRGV10*01"),
        cdr3a = c("CATWDRF",   "CATWDRF",   "CATWDSF"),
        vb    = c("TRDV1*01",  "TRDV2*01",  "TRDV1*01"),
        cdr3b = c("CALGELGDDKLIF", "CALGELGDDKLIF", "CALGELSDDKLIF"),
        stringsAsFactors = FALSE
    )

    # Human GD
    obj <- TCRrep(gd_df, "human_gd", chains = "GD")
    expect_true(is(obj, "TCRrep"))
    expect_equal(obj@chains, "GD")
    expect_equal(obj@organism, "human_gd")

    # Human GD with distances
    obj_dist <- TCRrep(gd_df, "human_gd", chains = "GD",
                        compute_distances = TRUE)
    expect_false(is.null(obj_dist@paired_dist))
    expect_equal(dim(obj_dist@paired_dist), c(3L, 3L))
    expect_equal(unname(diag(obj_dist@paired_dist)), c(0, 0, 0))

    # Mouse GD
    gd_mouse <- data.frame(
        va = c("TRGV1*01", "TRGV1*01"), cdr3a = c("CATWDRF", "CATWDSF"),
        vb = c("TRAV1*01", "TRAV1*01"),
        cdr3b = c("CALGELGDDKLIF", "CALGELSDDKLIF"),
        stringsAsFactors = FALSE
    )
    obj_mouse <- TCRrep(gd_mouse, "mouse_gd", chains = "GD")
    expect_equal(obj_mouse@organism, "mouse_gd")
})
