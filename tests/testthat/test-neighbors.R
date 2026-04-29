test_that("tcrdist_knn returns correctly shaped output", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    knn <- tcrdist_knn(tcrs, "human", K = 1L)
    expect_true(is.list(knn))
    expect_true(!is.null(knn$knn_indices))
    expect_true(!is.null(knn$knn_distances))
    expect_equal(dim(knn$knn_indices),   c(3L, 1L))
    expect_equal(dim(knn$knn_distances), c(3L, 1L))
})

test_that("tcrdist_knn indices are 1-based and within range", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    knn <- tcrdist_knn(tcrs, "human", K = 1L)
    idx <- as.vector(knn$knn_indices)
    expect_true(all(idx >= 1L))
    expect_true(all(idx <= 3L))
})

test_that("tcrdist_knn K=2 sorted by distance when sort_nbrs=TRUE", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    knn <- tcrdist_knn(tcrs, "human", K = 2L, sort_nbrs = TRUE)
    # Each row should have non-decreasing distances
    for (i in seq_len(nrow(knn$knn_distances))) {
        dists <- knn$knn_distances[i, ]
        expect_true(dists[1] <= dists[2])
    }
})

test_that("knn_from_matrix and tcrdist_knn agree on K=1 nearest neighbors", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    mat      <- tcrdist_matrix(tcrs, "human")
    knn_tcr  <- tcrdist_knn(tcrs, "human", K = 1L)
    knn_mat  <- knn_from_matrix(mat, K = 1L)
    expect_equal(knn_tcr$knn_indices,   knn_mat$knn_indices)
    expect_equal(knn_tcr$knn_distances, knn_mat$knn_distances)
})

test_that("tcrdist_radius_neighbors returns TCR 1 and 2 as neighbors at radius=50", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    # dist(1,2)=28 <= 50; dist(1,3)=239 > 50; dist(2,3)=227 > 50
    nbrs <- tcrdist_radius_neighbors(tcrs, "human", radius = 50)
    expect_equal(length(nbrs), 3L)
    # TCR 1: neighbor is TCR 2
    expect_true(2L %in% nbrs[[1]]$indices)
    # TCR 2: neighbor is TCR 1
    expect_true(1L %in% nbrs[[2]]$indices)
    # TCR 3: no neighbors within radius 50
    expect_equal(length(nbrs[[3]]$indices), 0L)
})

test_that("tcrdist_radius_neighbors all returned distances are within radius", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    radius <- 250
    nbrs   <- tcrdist_radius_neighbors(tcrs, "human", radius = radius)
    for (i in seq_along(nbrs)) {
        if (length(nbrs[[i]]$distances) > 0L) {
            expect_true(all(nbrs[[i]]$distances <= radius))
        }
    }
})

test_that("tcrdist_radius_neighbors returns 1-based indices", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    nbrs <- tcrdist_radius_neighbors(tcrs, "human", radius = 250)
    n <- nrow(tcrs)
    for (i in seq_along(nbrs)) {
        idx <- nbrs[[i]]$indices
        if (length(idx) > 0L) {
            expect_true(all(idx >= 1L))
            expect_true(all(idx <= n))
        }
    }
})

test_that("group masking with all same group produces sentinel distances", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV20-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )
    # All TCRs in the same group: every pair is masked -> sentinel 1e3
    knn <- tcrdist_knn(tcrs, "human", K = 1L,
                       agroups = c(1L, 1L, 1L),
                       bgroups = c(1L, 1L, 1L))
    expect_true(all(knn$knn_distances == 1e3))
})

test_that("knn_from_matrix rejects non-square matrix", {
    D <- matrix(1:6, nrow = 2L, ncol = 3L)
    expect_error(knn_from_matrix(D, K = 1L), "square")
})

test_that("knn_from_matrix rejects K >= N", {
    D <- matrix(c(0, 1, 1, 0), nrow = 2L)
    expect_error(knn_from_matrix(D, K = 2L), "K must be")
})

test_that("knn_from_pca returns correct output shape", {
    set.seed(42L)
    pca <- matrix(rnorm(20L), nrow = 5L, ncol = 4L)
    knn <- knn_from_pca(pca, K = 2L)
    expect_equal(dim(knn$knn_indices),   c(5L, 2L))
    expect_equal(dim(knn$knn_distances), c(5L, 2L))
})

test_that("knn_from_pca returns 1-based indices", {
    set.seed(42L)
    pca <- matrix(rnorm(20L), nrow = 5L, ncol = 4L)
    knn <- knn_from_pca(pca, K = 2L)
    idx <- as.vector(knn$knn_indices)
    expect_true(all(idx >= 1L))
    expect_true(all(idx <= 5L))
})

test_that("tcrdist_knn rejects K out of range", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
        stringsAsFactors = FALSE
    )
    expect_error(tcrdist_knn(tcrs, "human", K = 2L), "K must")
    expect_error(tcrdist_knn(tcrs, "human", K = 0L), "K must")
})

test_that("tcrdist_radius_neighbors rejects negative radius", {
    tcrs <- data.frame(
        va    = c("TRAV1-1*01"),
        cdr3a = c("CAVRDSSYKLIF"),
        vb    = c("TRBV19*01"),
        cdr3b = c("CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )
    expect_error(tcrdist_radius_neighbors(tcrs, "human", radius = -1), "non-negative")
})
