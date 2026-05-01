# Tests for single-chain (alpha-only / beta-only) TCRdist support.

# -- Shared test data --------------------------------------------------------

make_paired_tcrs <- function() {
    data.frame(
        va    = c("TRAV1-1*01", "TRAV1-2*01", "TRAV12-1*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CALSDRSYEKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV5-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CASSLGTEAFF"),
        stringsAsFactors = FALSE
    )
}


# -- .fill_missing_chain tests -----------------------------------------------

test_that(".fill_missing_chain returns AB for paired input", {
    tcrs <- make_paired_tcrs()
    result <- tcrdistR:::.fill_missing_chain(tcrs, "human")
    expect_equal(result$chain, "AB")
    expect_equal(ncol(result$tcrs), ncol(tcrs))
})

test_that(".fill_missing_chain fills alpha for beta-only input", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]
    result <- tcrdistR:::.fill_missing_chain(beta_only, "human")
    expect_equal(result$chain, "B")
    expect_true("va" %in% colnames(result$tcrs))
    expect_true("cdr3a" %in% colnames(result$tcrs))
    expect_true(all(!is.na(result$tcrs$va)))
    expect_true(all(!is.na(result$tcrs$cdr3a)))
})

test_that(".fill_missing_chain fills beta for alpha-only input", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]
    result <- tcrdistR:::.fill_missing_chain(alpha_only, "human")
    expect_equal(result$chain, "A")
    expect_true("vb" %in% colnames(result$tcrs))
    expect_true("cdr3b" %in% colnames(result$tcrs))
    expect_true(all(!is.na(result$tcrs$vb)))
    expect_true(all(!is.na(result$tcrs$cdr3b)))
})

test_that(".fill_missing_chain errors with no chain columns", {
    tcrs <- data.frame(x = 1:3)
    expect_error(
        tcrdistR:::.fill_missing_chain(tcrs, "human"),
        "at least alpha.*or beta"
    )
})

test_that(".fill_missing_chain preserves extra columns", {
    tcrs <- make_paired_tcrs()
    beta_only <- data.frame(
        vb = tcrs$vb, cdr3b = tcrs$cdr3b, epitope = c("A", "B", "C"),
        stringsAsFactors = FALSE
    )
    result <- tcrdistR:::.fill_missing_chain(beta_only, "human")
    expect_true("epitope" %in% colnames(result$tcrs))
    expect_equal(result$tcrs$epitope, c("A", "B", "C"))
})


# -- tcrdist_matrix single-chain tests --------------------------------------
# Component additivity (single-chain == paired+components) is tested in
# test-components.R.  Here we test that single-chain input produces a valid
# distance matrix.

test_that("beta-only tcrdist_matrix produces valid output", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]
    d <- tcrdist_matrix(beta_only, "human")
    expect_equal(nrow(d), 3L)
    expect_equal(ncol(d), 3L)
})

test_that("alpha-only tcrdist_matrix produces valid output", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]
    d <- tcrdist_matrix(alpha_only, "human")
    expect_equal(nrow(d), 3L)
    expect_equal(ncol(d), 3L)
})

# -- tcrdist_sparse single-chain tests --------------------------------------
# Sparse/rect additivity vs paired+components is tested in test-components.R.

test_that("beta-only tcrdist_sparse produces valid output", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]
    sp <- tcrdist_sparse(beta_only, "human", threshold = Inf)
    expect_true(inherits(sp, "dgCMatrix"))
    expect_equal(nrow(sp), 3L)
})

test_that("alpha-only tcrdist_sparse produces valid output", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]
    sp <- tcrdist_sparse(alpha_only, "human", threshold = Inf)
    expect_true(inherits(sp, "dgCMatrix"))
    expect_equal(nrow(sp), 3L)
})


# -- tcrdist_rect single-chain tests ----------------------------------------

test_that("beta-only tcrdist_rect produces correct dimensions", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]
    ref <- beta_only[c(1, 3), ]
    r <- tcrdist_rect(beta_only, ref, "human")
    expect_equal(dim(r), c(3L, 2L))
})

test_that("alpha-only tcrdist_rect produces correct dimensions", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]
    ref <- alpha_only[c(1, 3), ]
    r <- tcrdist_rect(alpha_only, ref, "human")
    expect_equal(dim(r), c(3L, 2L))
})

test_that("tcrdist_rect errors when query and ref have different chains", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]
    beta_only  <- tcrs[, c("vb", "cdr3b")]
    expect_error(
        tcrdist_rect(alpha_only, beta_only, "human"),
        "same chain columns"
    )
})


# -- tcrdist_knn / radius_neighbors single-chain tests ----------------------
# KNN-vs-matrix consistency is tested in test-neighbors.R.

test_that("beta-only tcrdist_knn produces valid output", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]
    knn <- tcrdist_knn(beta_only, "human", K = 2L)
    expect_equal(nrow(knn$knn_indices), 3L)
    expect_equal(ncol(knn$knn_indices), 2L)
})

test_that("alpha-only tcrdist_knn produces valid output", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]
    knn <- tcrdist_knn(alpha_only, "human", K = 1L)
    expect_equal(nrow(knn$knn_indices), 3L)
    expect_equal(ncol(knn$knn_indices), 1L)
})

test_that("beta-only tcrdist_radius_neighbors produces valid output", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]
    nbrs <- tcrdist_radius_neighbors(beta_only, "human", radius = 200)
    expect_length(nbrs, 3L)
    expect_true(all(vapply(nbrs, function(x) is.list(x), logical(1))))
})


# -- TCRrep single-chain tests ----------------------------------------------

test_that("TCRrep with chains='B' and compute_distances works", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]

    rep_b <- TCRrep(beta_only, organism = "human", chains = "B",
                    compute_distances = TRUE, deduplicate = FALSE)
    expect_false(is.null(rep_b@paired_dist))
    expect_equal(nrow(rep_b@paired_dist), 3L)
    expect_equal(ncol(rep_b@paired_dist), 3L)
    expect_true(all(diag(rep_b@paired_dist) == 0))
})

test_that("TCRrep with chains='A' and compute_distances works", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]

    rep_a <- TCRrep(alpha_only, organism = "human", chains = "A",
                    compute_distances = TRUE, deduplicate = FALSE)
    expect_false(is.null(rep_a@paired_dist))
    expect_equal(nrow(rep_a@paired_dist), 3L)
})

test_that("TCRrep single-chain distances match paired components", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]

    rep_b <- TCRrep(beta_only, organism = "human", chains = "B",
                    compute_distances = TRUE, deduplicate = FALSE)
    d_ref <- tcrdist_matrix(tcrs, "human", components = "beta")
    expect_equal(rep_b@paired_dist, d_ref, tolerance = 1e-10)
})


# -- Mouse organism tests ---------------------------------------------------

test_that("single-chain works with mouse organism", {
    tcrs <- data.frame(
        vb    = c("TRBV13-1*01", "TRBV13-1*01", "TRBV29*01"),
        cdr3b = c("CASSDAGGRNTLYF", "CASSDAGGNTLYF", "CASSPDRGEVFF"),
        stringsAsFactors = FALSE
    )
    d <- tcrdist_matrix(tcrs, "mouse")
    expect_equal(nrow(d), 3L)
    expect_equal(ncol(d), 3L)
    expect_true(all(diag(d) == 0))
})
