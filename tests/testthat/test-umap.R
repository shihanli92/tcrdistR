test_that("compute_tcrdist_umap works from tcr_df + organism", {
    skip_if_not_installed("uwot")
    data(dash)
    sub <- dash[1:100, ]

    result <- compute_tcrdist_umap(sub, organism = "mouse", seed = 42)
    expect_true(is.list(result))
    expect_true(is.matrix(result$embeddings))
    expect_equal(nrow(result$embeddings), 100L)
    expect_equal(ncol(result$embeddings), 2L)
    expect_equal(result$n_components, 2L)
    expect_true(is.matrix(result$pca_embeddings))
})

test_that("compute_tcrdist_umap works from pca_embeddings", {
    skip_if_not_installed("uwot")
    data(dash)
    sub <- dash[1:100, ]

    pca <- compute_tcrdist_kernel_pca(sub, "mouse", n_components = 20L)
    result <- compute_tcrdist_umap(pca_embeddings = pca$embeddings, seed = 42)

    expect_equal(nrow(result$embeddings), 100L)
    expect_equal(ncol(result$embeddings), 2L)
    expect_equal(nrow(result$pca_embeddings), 100L)
    expect_equal(ncol(result$pca_embeddings), 20L)
})

test_that("compute_tcrdist_umap n_components = 3 produces 3D output", {
    skip_if_not_installed("uwot")

    mat <- matrix(rnorm(500), nrow = 50, ncol = 10)
    result <- compute_tcrdist_umap(pca_embeddings = mat, n_components = 3L,
                                    seed = 1)
    expect_equal(ncol(result$embeddings), 3L)
    expect_equal(result$n_components, 3L)
})

test_that("compute_tcrdist_umap is reproducible with seed", {
    skip_if_not_installed("uwot")

    mat <- matrix(rnorm(300), nrow = 30, ncol = 10)
    r1 <- compute_tcrdist_umap(pca_embeddings = mat, seed = 123)
    r2 <- compute_tcrdist_umap(pca_embeddings = mat, seed = 123)
    expect_equal(r1$embeddings, r2$embeddings)
})

test_that("compute_tcrdist_umap errors without inputs", {
    skip_if_not_installed("uwot")

    expect_error(
        compute_tcrdist_umap(),
        "provide either"
    )
})

test_that("compute_tcrdist_umap errors when tcr_df given without organism", {
    skip_if_not_installed("uwot")
    data(dash)

    expect_error(
        compute_tcrdist_umap(tcr_df = dash[1:10, ]),
        "organism.*required"
    )
})

test_that("compute_tcrdist_umap clamps n_neighbors to dataset size", {
    skip_if_not_installed("uwot")

    # Small dataset: n_neighbors > n-1 should not error
    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    result <- compute_tcrdist_umap(pca_embeddings = mat, n_neighbors = 100L,
                                    seed = 1)
    expect_equal(nrow(result$embeddings), 5L)
})
