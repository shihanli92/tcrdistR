# Tests for hierarchical clustering and neighborhood tests


# ===========================================================================
# tcrdist_hclust
# ===========================================================================

test_that("tcrdist_hclust returns correct structure", {
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
    # hclust has n-1 merges
    expect_equal(nrow(result$hclust$merge), 4L)
})

test_that("tcrdist_hclust subsamples when n > max_tcrs", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = rep("TRAV1-1*01", 10),
        cdr3a = paste0("CAVRDSSYK", sample(c("A","C","D","E","F","G","H","I","K","L"), 10)),
        vb = rep("TRBV5-1*01", 10),
        cdr3b = paste0("CASSIRSSYE", sample(c("A","C","D","E","F","G","H","I","K","L"), 10)),
        stringsAsFactors = FALSE
    )

    result <- tcrdist_hclust(tcr_df, "human", max_tcrs = 5L)
    expect_equal(length(result$indices), 5L)
    expect_equal(nrow(result$dist_matrix), 5L)
})


# ===========================================================================
# cluster_tcrs
# ===========================================================================

test_that("cluster_tcrs with k gives correct number of clusters", {
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

    clusters <- cluster_tcrs(tcr_df, "human", k = 3)
    expect_length(clusters, 5L)
    expect_equal(length(unique(clusters)), 3L)
})

test_that("cluster_tcrs with h assigns all TCRs", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01", "TRAV10*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01", "TRBV7-2*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY"),
        stringsAsFactors = FALSE
    )

    clusters <- cluster_tcrs(tcr_df, "human", h = 50)
    expect_length(clusters, 3L)
    expect_true(all(clusters >= 1L))
})

test_that("cluster_tcrs errors without k or h", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY"),
        stringsAsFactors = FALSE
    )

    expect_error(cluster_tcrs(tcr_df, "human"), "exactly one")
})

test_that("cluster_tcrs with dist_matrix param works for hierarchical", {
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

    dm <- tcrdist_matrix(tcr_df, "human")
    cl1 <- cluster_tcrs(tcr_df, "human", k = 3)
    cl2 <- cluster_tcrs(dist_matrix = dm, k = 3)
    expect_equal(cl1, cl2)
})


# ===========================================================================
# cluster_tcrs: Leiden
# ===========================================================================

test_that("cluster_tcrs method='leiden' returns valid clusters", {
    skip_on_cran()
    skip_if_not_installed("igraph")

    data(dash, envir = environment())
    sub <- dash[1:30, ]

    cl <- cluster_tcrs(sub, "mouse", method = "leiden")
    expect_length(cl, 30L)
    expect_true(all(cl >= 1L))
    expect_true(is.integer(cl))
})


# ===========================================================================
# cluster_tcrs: Louvain
# ===========================================================================

test_that("cluster_tcrs method='louvain' returns valid clusters", {
    skip_on_cran()
    skip_if_not_installed("igraph")

    data(dash, envir = environment())
    sub <- dash[1:30, ]

    cl <- cluster_tcrs(sub, "mouse", method = "louvain")
    expect_length(cl, 30L)
    expect_true(all(cl >= 1L))
    expect_true(is.integer(cl))
})

test_that("cluster_tcrs leiden/louvain requires tcr_df and organism", {
    skip_on_cran()
    skip_if_not_installed("igraph")

    dm <- matrix(0, nrow = 5, ncol = 5)
    expect_error(cluster_tcrs(dist_matrix = dm, method = "leiden"),
                 "tcr_df.*organism.*required")
})


# ===========================================================================
# cluster_tcrs: DBSCAN
# ===========================================================================

test_that("cluster_tcrs method='dbscan' with explicit eps", {
    skip_on_cran()
    skip_if_not_installed("dbscan")

    data(dash, envir = environment())
    sub <- dash[1:30, ]

    cl <- cluster_tcrs(sub, "mouse", method = "dbscan", eps = 50)
    expect_length(cl, 30L)
    expect_true(all(cl >= 0L))
    expect_true(is.integer(cl))
})

test_that("cluster_tcrs method='dbscan' with auto eps", {
    skip_on_cran()
    skip_if_not_installed("dbscan")

    data(dash, envir = environment())
    sub <- dash[1:30, ]

    expect_message(
        cl <- cluster_tcrs(sub, "mouse", method = "dbscan"),
        "Auto-detected DBSCAN eps"
    )
    expect_length(cl, 30L)
    expect_true(all(cl >= 0L))
})

test_that("cluster_tcrs dbscan works with precomputed dist_matrix", {
    skip_on_cran()
    skip_if_not_installed("dbscan")

    data(dash, envir = environment())
    sub <- dash[1:30, ]
    dm <- tcrdist_matrix(sub, "mouse")

    cl <- cluster_tcrs(dist_matrix = dm, method = "dbscan", eps = 50)
    expect_length(cl, 30L)
    expect_true(all(cl >= 0L))
})


# ===========================================================================
# cluster_tcrs: K-medoids
# ===========================================================================

test_that("cluster_tcrs method='kmedoids' returns correct k clusters", {
    skip_on_cran()
    skip_if_not_installed("cluster")

    data(dash, envir = environment())
    sub <- dash[1:30, ]

    cl <- cluster_tcrs(sub, "mouse", method = "kmedoids", k = 4)
    expect_length(cl, 30L)
    expect_equal(length(unique(cl)), 4L)
    expect_true(all(cl >= 1L))
    expect_true(is.integer(cl))
})

test_that("cluster_tcrs kmedoids errors without k", {
    skip_on_cran()
    skip_if_not_installed("cluster")

    data(dash, envir = environment())
    sub <- dash[1:10, ]

    expect_error(cluster_tcrs(sub, "mouse", method = "kmedoids"),
                 "'k' is required")
})

test_that("cluster_tcrs kmedoids works with precomputed dist_matrix", {
    skip_on_cran()
    skip_if_not_installed("cluster")

    data(dash, envir = environment())
    sub <- dash[1:30, ]
    dm <- tcrdist_matrix(sub, "mouse")

    cl1 <- cluster_tcrs(sub, "mouse", method = "kmedoids", k = 3)
    cl2 <- cluster_tcrs(dist_matrix = dm, method = "kmedoids", k = 3)
    expect_equal(cl1, cl2)
})


# ===========================================================================
# .neighborhood_tally
# ===========================================================================

test_that(".neighborhood_tally: row sums match neighborhood sizes", {
    dist_mat <- matrix(c(0, 5, 100,
                          5, 0, 8,
                          100, 8, 0), nrow = 3)
    membership <- factor(c("A", "B", "B"))

    tally <- tcrdistR:::.neighborhood_tally(dist_mat, membership, radius = 10)
    expect_equal(nrow(tally), 3L)
    expect_equal(ncol(tally), 2L)  # A and B
    # Row 1: neighbors within 10 are {1, 2} -> 1 A + 1 B = 2
    expect_equal(sum(tally[1, ]), 2L)
    # Row 2: neighbors within 10 are {1, 2, 3} -> 1 A + 2 B = 3
    expect_equal(sum(tally[2, ]), 3L)
    # Row 3: neighbors within 10 are {2, 3} -> 0 A + 2 B = 2
    expect_equal(sum(tally[3, ]), 2L)
})


# ===========================================================================
# neighborhood_test
# ===========================================================================

test_that("neighborhood_test returns valid p-values", {
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

    result <- neighborhood_test(tcr_df, "human", variable,
                                 radius = 100, test = "fisher")
    expect_s3_class(result, "data.frame")
    expect_true(all(result$p_value >= 0 & result$p_value <= 1))
    expect_true(all(result$p_adjusted >= 0 & result$p_adjusted <= 1))
    expect_true(all(result$odds_ratio >= 0 | is.na(result$odds_ratio)))
    expect_equal(nrow(result), 5L)
})

test_that("neighborhood_test: fisher requires binary variable", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01", "TRAV10*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01", "TRBV7-2*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY"),
        stringsAsFactors = FALSE
    )

    expect_error(
        neighborhood_test(tcr_df, "human", c("A", "B", "C"),
                           test = "fisher"),
        "exactly 2"
    )
})


# ===========================================================================
# dist_matrix bypass
# ===========================================================================

test_that("tcrdist_hclust accepts precomputed dist_matrix", {
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

    dm <- tcrdist_matrix(tcr_df, "human")
    result <- tcrdist_hclust(dist_matrix = dm)
    expect_s3_class(result$hclust, "hclust")
    expect_equal(nrow(result$dist_matrix), 5L)
})

test_that("tcrdist_hclust dist_matrix matches tcr_df path", {
    skip_on_cran()

    tcr_df <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-2*01", "TRAV10*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF"),
        vb = c("TRBV5-1*01", "TRBV6-1*01", "TRBV7-2*01"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY"),
        stringsAsFactors = FALSE
    )

    dm <- tcrdist_matrix(tcr_df, "human")
    r1 <- tcrdist_hclust(tcr_df, "human")
    r2 <- tcrdist_hclust(dist_matrix = dm)
    expect_equal(r1$dist_matrix, r2$dist_matrix)
})

test_that("neighborhood_test accepts precomputed dist_matrix", {
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

    dm <- tcrdist_matrix(tcr_df, "human")
    variable <- c("A", "B", "A", "B", "A")
    r1 <- neighborhood_test(tcr_df, "human", variable, radius = 100,
                             test = "fisher")
    r2 <- neighborhood_test(variable = variable, radius = 100,
                             test = "fisher", dist_matrix = dm)
    expect_equal(r1$p_value, r2$p_value)
    expect_equal(r1$n_neighbors, r2$n_neighbors)
})
