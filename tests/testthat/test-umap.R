# ===========================================================================
# KNN path tests (tcr_df + organism)
# ===========================================================================

test_that("KNN path: basic output structure", {
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
    expect_equal(nrow(result$knn_distances), 100L)
    expect_true(inherits(result$knn_graph, "Matrix"))
    expect_equal(nrow(result$knn_graph), 100L)
    expect_equal(ncol(result$knn_graph), 100L)

    # Weighted NN distances
    expect_true(is.numeric(result$nndists))
    expect_equal(length(result$nndists), 100L)
    expect_true(all(result$nndists >= 0))

    # No clusters by default
    expect_null(result$clusters)
})

test_that("KNN path: n_components = 3 produces 3D output", {
    skip_if_not_installed("uwot")
    data(dash)
    sub <- dash[1:80, ]

    result <- compute_tcrdist_umap(sub, "mouse", n_components = 3L, seed = 1)
    expect_equal(ncol(result$embeddings), 3L)
    expect_equal(result$n_components, 3L)
})

test_that("KNN path: custom n_neighbors", {
    skip_if_not_installed("uwot")
    data(dash)
    sub <- dash[1:60, ]

    result <- compute_tcrdist_umap(sub, "mouse", n_neighbors = 10L, seed = 1)
    expect_equal(ncol(result$knn_indices), 10L)
    expect_equal(result$n_neighbors, 10L)
})

test_that("KNN path: reproducible with seed", {
    skip_if_not_installed("uwot")
    data(dash)
    sub <- dash[1:60, ]

    r1 <- compute_tcrdist_umap(sub, "mouse", seed = 123)
    r2 <- compute_tcrdist_umap(sub, "mouse", seed = 123)
    expect_equal(r1$embeddings, r2$embeddings)
    expect_equal(r1$knn_indices, r2$knn_indices)
    expect_equal(r1$knn_distances, r2$knn_distances)
})

test_that("KNN path: clustering produces integer vector", {
    skip_if_not_installed("uwot")
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:100, ]

    result <- compute_tcrdist_umap(sub, "mouse", seed = 42, cluster = TRUE)
    expect_true(is.integer(result$clusters))
    expect_equal(length(result$clusters), 100L)
    expect_true(all(result$clusters >= 0L))
})

test_that("KNN path: knn_graph is symmetric sparse matrix", {
    skip_if_not_installed("uwot")
    data(dash)
    sub <- dash[1:50, ]

    result <- compute_tcrdist_umap(sub, "mouse", seed = 1)
    g <- result$knn_graph
    expect_true(inherits(g, "dgCMatrix") || inherits(g, "Matrix"))
    expect_true(Matrix::isSymmetric(g, tol = 1e-10))
})

# ===========================================================================
# PCA path tests (pca_embeddings)
# ===========================================================================

test_that("PCA path: works from pca_embeddings", {
    skip_if_not_installed("uwot")
    data(dash)
    sub <- dash[1:100, ]

    pca <- compute_tcrdist_kernel_pca(sub, "mouse", n_components = 20L)
    result <- compute_tcrdist_umap(pca_embeddings = pca$embeddings, seed = 42)

    expect_equal(result$method, "pca")
    expect_equal(nrow(result$embeddings), 100L)
    expect_equal(ncol(result$embeddings), 2L)
    expect_equal(nrow(result$pca_embeddings), 100L)
    expect_equal(ncol(result$pca_embeddings), 20L)
})

test_that("PCA path: n_components = 3 produces 3D output", {
    skip_if_not_installed("uwot")

    mat <- matrix(rnorm(500), nrow = 50, ncol = 10)
    result <- compute_tcrdist_umap(pca_embeddings = mat, n_components = 3L,
                                    seed = 1)
    expect_equal(ncol(result$embeddings), 3L)
    expect_equal(result$n_components, 3L)
})

test_that("PCA path: reproducible with seed", {
    skip_if_not_installed("uwot")

    mat <- matrix(rnorm(300), nrow = 30, ncol = 10)
    r1 <- compute_tcrdist_umap(pca_embeddings = mat, seed = 123)
    r2 <- compute_tcrdist_umap(pca_embeddings = mat, seed = 123)
    expect_equal(r1$embeddings, r2$embeddings)
})

test_that("PCA path: clamps n_neighbors to dataset size", {
    skip_if_not_installed("uwot")

    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    result <- compute_tcrdist_umap(pca_embeddings = mat, n_neighbors = 100L,
                                    seed = 1)
    expect_equal(nrow(result$embeddings), 5L)
})

# ===========================================================================
# Error handling
# ===========================================================================

test_that("errors without inputs", {
    skip_if_not_installed("uwot")
    expect_error(compute_tcrdist_umap(), "provide either")
})

test_that("errors when tcr_df given without organism", {
    skip_if_not_installed("uwot")
    data(dash)
    expect_error(
        compute_tcrdist_umap(tcr_df = dash[1:10, ]),
        "organism.*required"
    )
})

# ===========================================================================
# Internal helper tests
# ===========================================================================

test_that(".smooth_knn_dist returns rho and sigma of correct length", {
    dmat <- matrix(c(1, 2, 3, 4, 0.5, 1.5, 2.5, 3.5), nrow = 2, byrow = TRUE)
    result <- tcrdistR:::.smooth_knn_dist(dmat)
    expect_equal(length(result$rho), 2L)
    expect_equal(length(result$sigma), 2L)
    expect_true(all(result$sigma > 0))
})

test_that(".calc_nndists returns correct weighted distances", {
    # K=3, weights: 1.0, 0.667, 0.333 (normalized)
    dmat <- matrix(c(1, 2, 3), nrow = 1)
    result <- tcrdistR:::.calc_nndists(dmat)
    expect_equal(length(result), 1L)
    expect_true(result > 0)

    # Empty input
    expect_equal(length(tcrdistR:::.calc_nndists(matrix(0, nrow = 0, ncol = 3))), 0L)
})

test_that(".build_knn_graph produces symmetric sparse matrix", {
    # Simple 4-point KNN with K=2
    idx <- matrix(c(2L, 3L, 1L, 3L, 1L, 4L, 3L, 2L), nrow = 4, byrow = TRUE)
    dst <- matrix(c(1.0, 2.0, 1.0, 1.5, 2.0, 0.5, 0.5, 1.5), nrow = 4, byrow = TRUE)

    g <- tcrdistR:::.build_knn_graph(idx, dst, 4L)
    expect_true(inherits(g, "Matrix"))
    expect_equal(nrow(g), 4L)
    expect_equal(ncol(g), 4L)
    expect_true(Matrix::isSymmetric(g, tol = 1e-10))
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

    # Should complete without infinite loop even with small n_neighbors
    result <- expect_no_error(
        compute_tcrdist_umap(tcr_df, "human", n_neighbors = 2L, seed = 42L)
    )
    expect_equal(nrow(result$embeddings), 6L)
})
