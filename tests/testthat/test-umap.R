# Tests for UMAP dimensionality reduction.
# Internal helpers (.smooth_knn_dist, .calc_nndists, .build_knn_graph)
# are tested implicitly through compute_tcrdist_umap.


# ===========================================================================
# KNN path (tcr_df + organism)
# ===========================================================================

test_that("KNN path: output structure, KNN artifacts, and no default clusters", {
    skip_if_not_installed("uwot")
    data(dash)
    sub <- dash[1:100, ]

    result <- compute_tcrdist_umap(sub, organism = "mouse", seed = 42)

    expect_true(is.list(result))
    expect_equal(result$method, "knn")
    expect_true(is.matrix(result$embeddings))
    expect_equal(nrow(result$embeddings), 100L)
    expect_equal(ncol(result$embeddings), 2L)
    expect_equal(result$n_components, 2L)
    # KNN artifacts
    expect_true(is.matrix(result$knn_indices))
    expect_true(is.matrix(result$knn_distances))
    expect_equal(nrow(result$knn_indices), 100L)
    expect_true(inherits(result$knn_graph, "Matrix"))
    expect_true(Matrix::isSymmetric(result$knn_graph, tol = 1e-10))
    # Weighted NN distances
    expect_equal(length(result$nndists), 100L)
    expect_true(all(result$nndists >= 0))
    # No clusters by default
    expect_null(result$clusters)
})

test_that("KNN path: n_components, n_neighbors, seed, and clustering", {
    skip_if_not_installed("uwot")
    data(dash)
    sub <- dash[1:80, ]

    # 3D
    r3d <- compute_tcrdist_umap(sub, "mouse", n_components = 3L, seed = 1)
    expect_equal(ncol(r3d$embeddings), 3L)

    # Custom n_neighbors
    r_nn <- compute_tcrdist_umap(dash[1:60, ], "mouse",
                                  n_neighbors = 10L, seed = 1)
    expect_equal(ncol(r_nn$knn_indices), 10L)

    # Reproducibility
    r1 <- compute_tcrdist_umap(dash[1:60, ], "mouse", seed = 123)
    r2 <- compute_tcrdist_umap(dash[1:60, ], "mouse", seed = 123)
    expect_equal(r1$embeddings, r2$embeddings)

    # Clustering
    skip_if_not_installed("igraph")
    r_cl <- compute_tcrdist_umap(dash[1:100, ], "mouse",
                                  seed = 42, cluster = TRUE)
    expect_true(is.integer(r_cl$clusters))
    expect_equal(length(r_cl$clusters), 100L)
    expect_true(all(r_cl$clusters >= 0L))
})


# ===========================================================================
# PCA path (pca_embeddings)
# ===========================================================================

test_that("PCA path: works with various options", {
    skip_if_not_installed("uwot")
    data(dash)

    pca <- compute_tcrdist_kernel_pca(dash[1:100, ], "mouse",
                                       n_components = 20L)
    result <- compute_tcrdist_umap(pca_embeddings = pca$embeddings, seed = 42)
    expect_equal(result$method, "pca")
    expect_equal(nrow(result$embeddings), 100L)
    expect_equal(ncol(result$embeddings), 2L)

    # 3D
    mat <- matrix(rnorm(500), nrow = 50, ncol = 10)
    r3d <- compute_tcrdist_umap(pca_embeddings = mat, n_components = 3L,
                                 seed = 1)
    expect_equal(ncol(r3d$embeddings), 3L)

    # Reproducibility
    mat2 <- matrix(rnorm(300), nrow = 30, ncol = 10)
    r1 <- compute_tcrdist_umap(pca_embeddings = mat2, seed = 123)
    r2 <- compute_tcrdist_umap(pca_embeddings = mat2, seed = 123)
    expect_equal(r1$embeddings, r2$embeddings)

    # Clamps n_neighbors to dataset size
    mat3 <- matrix(rnorm(50), nrow = 5, ncol = 10)
    r_clamp <- compute_tcrdist_umap(pca_embeddings = mat3,
                                     n_neighbors = 100L, seed = 1)
    expect_equal(nrow(r_clamp$embeddings), 5L)
})


# ===========================================================================
# Error handling and edge cases
# ===========================================================================

test_that("compute_tcrdist_umap rejects bad inputs", {
    skip_if_not_installed("uwot")
    expect_error(compute_tcrdist_umap(), "provide either")
    data(dash)
    expect_error(compute_tcrdist_umap(tcr_df = dash[1:10, ]),
                 "organism.*required")
})

test_that("compute_tcrdist_umap connectivity loop has bounded n_neighbors", {
    skip_on_cran()
    skip_if_not_installed("uwot")
    skip_if_not_installed("igraph")

    # Two very different TCR clusters — likely fragmented at small n_neighbors
    tcr_df <- data.frame(
        va = c(rep("TRAV1-1*01", 3), rep("TRAV12-2*01", 3)),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVRDSSYKLIF",
                   "CAVSANSGTYF", "CAVSANSGTYF", "CAVSANSGTYF"),
        vb = c(rep("TRBV19*01", 3), rep("TRBV20-1*01", 3)),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CASSIRSSYEQYF",
                   "CSARDRTGNTIYF", "CSARDRTGNTIYF", "CSARDRTGNTIYF"),
        stringsAsFactors = FALSE
    )

    result <- expect_no_error(
        compute_tcrdist_umap(tcr_df, "human", n_neighbors = 2L, seed = 42L)
    )
    expect_equal(nrow(result$embeddings), 6L)
})
