# Tests for hierarchical clustering and neighborhood tests


# ===========================================================================
# tcrdist_hclust
# ===========================================================================

test_that("tcrdist_hclust returns correct structure and supports subsampling", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01", "TRAV10*01",
               "TRAV1-1*01", "TRAV1-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF",
                   "CAVRDSSYKLIF", "CAVKDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01", "TRBV7-2*01",
               "TRBV5-1*01", "TRBV6-1*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY",
                   "CASSIRSSYEQF", "CASSIKSSYEQF"),
        stringsAsFactors = FALSE
    )

    result <- tcrdist_hclust(tcr_df, "human")
    expect_type(result, "list")
    expect_s3_class(result$hclust, "hclust")
    expect_true(is.matrix(result$dist_matrix))
    expect_equal(length(result$indices), 5L)
    expect_equal(nrow(result$hclust$merge), 4L)

    # Subsampling
    big_df <- data.frame(
        va = rep("TRAV1-1*01", 10),
        cdr3a = paste0("CAVRDSSYK", sample(c("A","C","D","E","F","G","H","I","K","L"), 10)),
        vb = rep("TRBV5-1*01", 10),
        cdr3b = paste0("CASSIRSSYE", sample(c("A","C","D","E","F","G","H","I","K","L"), 10)),
        stringsAsFactors = FALSE
    )
    result_sub <- tcrdist_hclust(big_df, "human", max_tcrs = 5L)
    expect_equal(length(result_sub$indices), 5L)
    expect_equal(nrow(result_sub$dist_matrix), 5L)

    # dist_matrix bypass matches tcr_df path
    dm <- tcrdist_matrix(tcr_df, "human")
    r_dm <- tcrdist_hclust(dist_matrix = dm)
    expect_s3_class(r_dm$hclust, "hclust")
    expect_equal(result$dist_matrix, r_dm$dist_matrix)
})


# ===========================================================================
# cluster_tcrs: hierarchical
# ===========================================================================

test_that("cluster_tcrs hierarchical with k, h, and dist_matrix bypass", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01", "TRAV10*01",
               "TRAV1-1*01", "TRAV1-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF",
                   "CAVRDSSYKLIF", "CAVKDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01", "TRBV7-2*01",
               "TRBV5-1*01", "TRBV6-1*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY",
                   "CASSIRSSYEQF", "CASSIKSSYEQF"),
        stringsAsFactors = FALSE
    )

    # k mode
    cl_k <- cluster_tcrs(tcr_df, "human", k = 3)
    expect_length(cl_k, 5L)
    expect_equal(length(unique(cl_k)), 3L)

    # h mode
    cl_h <- cluster_tcrs(tcr_df[1:3, ], "human", h = 50)
    expect_length(cl_h, 3L)
    expect_true(all(cl_h >= 1L))

    # Error without k or h
    expect_error(cluster_tcrs(tcr_df[1:2, ], "human"), "exactly one")

    # dist_matrix bypass matches tcr_df path
    dm <- tcrdist_matrix(tcr_df, "human")
    cl_dm <- cluster_tcrs(dist_matrix = dm, k = 3)
    expect_equal(cl_k, cl_dm)
})


# ===========================================================================
# cluster_tcrs: Leiden and Louvain
# ===========================================================================

test_that("cluster_tcrs leiden and louvain produce valid clusters", {
    skip_on_cran()
    skip_if_not_installed("igraph")

    data(dash, envir = environment())
    sub <- dash[1:30, ]

    cl_lei <- cluster_tcrs(sub, "mouse", method = "leiden")
    expect_length(cl_lei, 30L)
    expect_true(all(cl_lei >= 1L))
    expect_true(is.integer(cl_lei))

    cl_louv <- cluster_tcrs(sub, "mouse", method = "louvain")
    expect_length(cl_louv, 30L)
    expect_true(all(cl_louv >= 1L))
    expect_true(is.integer(cl_louv))

    # Requires tcr_df and organism (not dist_matrix only)
    dm <- matrix(0, nrow = 5, ncol = 5)
    expect_error(cluster_tcrs(dist_matrix = dm, method = "leiden"),
                 "tcr_df.*organism.*required")
})


# ===========================================================================
# cluster_tcrs: DBSCAN
# ===========================================================================

test_that("cluster_tcrs DBSCAN with explicit eps, auto eps, and dist_matrix", {
    skip_on_cran()
    skip_if_not_installed("dbscan")

    data(dash, envir = environment())
    sub <- dash[1:30, ]

    # Explicit eps
    cl_eps <- cluster_tcrs(sub, "mouse", method = "dbscan", eps = 50)
    expect_length(cl_eps, 30L)
    expect_true(all(cl_eps >= 0L))
    expect_true(is.integer(cl_eps))

    # Auto eps
    expect_message(
        cl_auto <- cluster_tcrs(sub, "mouse", method = "dbscan"),
        "Auto-detected DBSCAN eps"
    )
    expect_length(cl_auto, 30L)

    # dist_matrix bypass
    dm <- tcrdist_matrix(sub, "mouse")
    cl_dm <- cluster_tcrs(dist_matrix = dm, method = "dbscan", eps = 50)
    expect_length(cl_dm, 30L)
})


# ===========================================================================
# cluster_tcrs: K-medoids
# ===========================================================================

test_that("cluster_tcrs kmedoids with k and dist_matrix bypass", {
    skip_on_cran()
    skip_if_not_installed("cluster")

    data(dash, envir = environment())
    sub <- dash[1:30, ]

    cl <- cluster_tcrs(sub, "mouse", method = "kmedoids", k = 4)
    expect_length(cl, 30L)
    expect_equal(length(unique(cl)), 4L)
    expect_true(is.integer(cl))

    # Error without k
    expect_error(cluster_tcrs(sub, "mouse", method = "kmedoids"),
                 "'k' is required")

    # dist_matrix bypass
    dm <- tcrdist_matrix(sub, "mouse")
    cl_dm <- cluster_tcrs(dist_matrix = dm, method = "kmedoids", k = 4)
    expect_equal(cl, cl_dm)
})


# ===========================================================================
# neighborhood_test
# ===========================================================================

test_that("neighborhood_test returns valid results with dist_matrix bypass", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-1*01", "TRAV1-2*01",
               "TRAV1-2*01", "TRAV10*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVKDSSYKLIF",
                   "CAVKDSYKLIF", "CAVRDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV5-1*01", "TRBV6-1*01",
               "TRBV6-1*01", "TRBV7-2*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIRSSYEQF", "CASSIKSSYEQY",
                   "CASSIKSSYEQF", "CASSIRSYEQY"),
        stringsAsFactors = FALSE
    )

    variable <- c("X", "X", "Y", "Y", "X")

    r1 <- neighborhood_test(tcr_df, "human", variable,
                             radius = 100, test = "fisher")
    expect_s3_class(r1, "data.frame")
    expect_true(all(r1$p_value >= 0 & r1$p_value <= 1))
    expect_true(all(r1$p_adjusted >= 0 & r1$p_adjusted <= 1))
    expect_equal(nrow(r1), 5L)

    # dist_matrix bypass matches
    dm <- tcrdist_matrix(tcr_df, "human")
    r2 <- neighborhood_test(variable = variable, radius = 100,
                             test = "fisher", dist_matrix = dm)
    expect_equal(r1$p_value, r2$p_value)
    expect_equal(r1$n_neighbors, r2$n_neighbors)

    # Fisher requires binary variable
    expect_error(
        neighborhood_test(tcr_df[1:3, ], "human", c("A", "B", "C"),
                           test = "fisher"),
        "exactly 2"
    )
})
