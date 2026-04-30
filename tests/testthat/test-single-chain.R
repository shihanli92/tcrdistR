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

test_that("beta-only tcrdist_matrix matches components='beta' on paired data", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]

    d_single <- tcrdist_matrix(beta_only, "human")
    d_ref    <- tcrdist_matrix(tcrs, "human", components = "beta")
    expect_equal(d_single, d_ref, tolerance = 1e-10)
})

test_that("alpha-only tcrdist_matrix matches components='alpha' on paired data", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]

    d_single <- tcrdist_matrix(alpha_only, "human")
    d_ref    <- tcrdist_matrix(tcrs, "human", components = "alpha")
    expect_equal(d_single, d_ref, tolerance = 1e-10)
})

test_that("beta-only with explicit components='cdr3b' works", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]

    d_cdr3b <- tcrdist_matrix(beta_only, "human", components = "cdr3b")
    d_ref   <- tcrdist_matrix(tcrs, "human", components = "cdr3b")
    expect_equal(d_cdr3b, d_ref, tolerance = 1e-10)
})

test_that("alpha-only with explicit components='va' works", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]

    d_va  <- tcrdist_matrix(alpha_only, "human", components = "va")
    d_ref <- tcrdist_matrix(tcrs, "human", components = "va")
    expect_equal(d_va, d_ref, tolerance = 1e-10)
})

test_that("single-chain diagonal is zero", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]
    d <- tcrdist_matrix(beta_only, "human")
    expect_true(all(diag(d) == 0))
})

test_that("single-chain matrix is symmetric", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]
    d <- tcrdist_matrix(alpha_only, "human")
    expect_equal(d, t(d))
})


# -- tcrdist_sparse single-chain tests --------------------------------------

test_that("beta-only tcrdist_sparse matches paired components='beta'", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]

    sp_single <- as.matrix(tcrdist_sparse(beta_only, "human", threshold = Inf))
    sp_ref    <- as.matrix(tcrdist_sparse(tcrs, "human", threshold = Inf,
                                           components = "beta"))
    expect_equal(sp_single, sp_ref, tolerance = 1e-10)
})

test_that("alpha-only tcrdist_sparse matches paired components='alpha'", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]

    sp_single <- as.matrix(tcrdist_sparse(alpha_only, "human", threshold = Inf))
    sp_ref    <- as.matrix(tcrdist_sparse(tcrs, "human", threshold = Inf,
                                           components = "alpha"))
    expect_equal(sp_single, sp_ref, tolerance = 1e-10)
})


# -- tcrdist_rect single-chain tests ----------------------------------------

test_that("beta-only tcrdist_rect matches paired components='beta'", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]
    ref <- beta_only[c(1, 3), ]

    r_single <- tcrdist_rect(beta_only, ref, "human")
    r_ref    <- tcrdist_rect(tcrs, tcrs[c(1, 3), ], "human",
                              components = "beta")
    expect_equal(r_single, r_ref, tolerance = 1e-10)
})

test_that("alpha-only tcrdist_rect matches paired components='alpha'", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]
    ref <- alpha_only[c(1, 3), ]

    r_single <- tcrdist_rect(alpha_only, ref, "human")
    r_ref    <- tcrdist_rect(tcrs, tcrs[c(1, 3), ], "human",
                              components = "alpha")
    expect_equal(r_single, r_ref, tolerance = 1e-10)
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


# -- tcrdist_knn single-chain tests -----------------------------------------

test_that("beta-only tcrdist_knn consistent with tcrdist_matrix", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]

    d <- tcrdist_matrix(beta_only, "human")
    knn <- tcrdist_knn(beta_only, "human", K = 2L)

    for (i in seq_len(nrow(beta_only))) {
        for (k in seq_len(2L)) {
            j <- knn$knn_indices[i, k]
            expect_equal(knn$knn_distances[i, k], d[i, j], tolerance = 1e-10)
        }
    }
})

test_that("alpha-only tcrdist_knn works", {
    tcrs <- make_paired_tcrs()
    alpha_only <- tcrs[, c("va", "cdr3a")]
    knn <- tcrdist_knn(alpha_only, "human", K = 1L)
    expect_equal(nrow(knn$knn_indices), 3L)
    expect_equal(ncol(knn$knn_indices), 1L)
})


# -- tcrdist_radius_neighbors single-chain tests ----------------------------

test_that("beta-only tcrdist_radius_neighbors consistent with tcrdist_matrix", {
    tcrs <- make_paired_tcrs()
    beta_only <- tcrs[, c("vb", "cdr3b")]

    d <- tcrdist_matrix(beta_only, "human")
    nbrs <- tcrdist_radius_neighbors(beta_only, "human", radius = 200)

    for (i in seq_len(nrow(beta_only))) {
        idx <- nbrs[[i]]$indices
        dst <- nbrs[[i]]$distances
        for (k in seq_along(idx)) {
            expect_equal(dst[k], d[i, idx[k]], tolerance = 1e-10)
        }
    }
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
