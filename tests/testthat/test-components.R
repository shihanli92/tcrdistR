# Tests for per-component TCRdist distance selection (components parameter).

# -- Shared test data --------------------------------------------------------

make_test_tcrs <- function() {
    data.frame(
        va    = c("TRAV1-1*01", "TRAV1-2*01", "TRAV12-1*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CALSDRSYEKLIF"),
        vb    = c("TRBV19*01", "TRBV19*01", "TRBV5-1*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF", "CASSLGTEAFF"),
        stringsAsFactors = FALSE
    )
}


# -- .resolve_components tests -----------------------------------------------

test_that(".resolve_components returns all TRUE for 'all'", {
    comp <- tcrdistR:::.resolve_components("all")
    expect_true(comp$va)
    expect_true(comp$cdr3a)
    expect_true(comp$vb)
    expect_true(comp$cdr3b)
})

test_that(".resolve_components handles presets correctly", {
    cdr3 <- tcrdistR:::.resolve_components("cdr3")
    expect_false(cdr3$va)
    expect_true(cdr3$cdr3a)
    expect_false(cdr3$vb)
    expect_true(cdr3$cdr3b)

    v_reg <- tcrdistR:::.resolve_components("v_region")
    expect_true(v_reg$va)
    expect_false(v_reg$cdr3a)
    expect_true(v_reg$vb)
    expect_false(v_reg$cdr3b)

    alpha <- tcrdistR:::.resolve_components("alpha")
    expect_true(alpha$va)
    expect_true(alpha$cdr3a)
    expect_false(alpha$vb)
    expect_false(alpha$cdr3b)

    beta <- tcrdistR:::.resolve_components("beta")
    expect_false(beta$va)
    expect_false(beta$cdr3a)
    expect_true(beta$vb)
    expect_true(beta$cdr3b)
})

test_that(".resolve_components handles custom vector", {
    comp <- tcrdistR:::.resolve_components(c("cdr3a", "vb"))
    expect_false(comp$va)
    expect_true(comp$cdr3a)
    expect_true(comp$vb)
    expect_false(comp$cdr3b)
})

test_that(".resolve_components rejects unknown terms", {
    expect_error(
        tcrdistR:::.resolve_components("cdr1"),
        "Unknown components"
    )
    expect_error(
        tcrdistR:::.resolve_components(c("va", "xyz")),
        "Unknown components"
    )
})


# -- tcrdist_matrix additivity tests ----------------------------------------

test_that("v_region + cdr3 = all (tcrdist_matrix)", {
    tcrs  <- make_test_tcrs()
    d_all  <- tcrdist_matrix(tcrs, "human")
    d_v    <- tcrdist_matrix(tcrs, "human", components = "v_region")
    d_cdr3 <- tcrdist_matrix(tcrs, "human", components = "cdr3")
    expect_equal(d_v + d_cdr3, d_all, tolerance = 1e-10)
})

test_that("alpha + beta = all (tcrdist_matrix)", {
    tcrs <- make_test_tcrs()
    d_all <- tcrdist_matrix(tcrs, "human")
    d_a   <- tcrdist_matrix(tcrs, "human", components = "alpha")
    d_b   <- tcrdist_matrix(tcrs, "human", components = "beta")
    expect_equal(d_a + d_b, d_all, tolerance = 1e-10)
})

test_that("va + cdr3a + vb + cdr3b = all (tcrdist_matrix)", {
    tcrs   <- make_test_tcrs()
    d_all  <- tcrdist_matrix(tcrs, "human")
    d_va   <- tcrdist_matrix(tcrs, "human", components = "va")
    d_cdr3a <- tcrdist_matrix(tcrs, "human", components = "cdr3a")
    d_vb   <- tcrdist_matrix(tcrs, "human", components = "vb")
    d_cdr3b <- tcrdist_matrix(tcrs, "human", components = "cdr3b")
    expect_equal(d_va + d_cdr3a + d_vb + d_cdr3b, d_all, tolerance = 1e-10)
})

test_that("custom combo c(cdr3a, vb) = cdr3a + vb (tcrdist_matrix)", {
    tcrs    <- make_test_tcrs()
    d_cdr3a <- tcrdist_matrix(tcrs, "human", components = "cdr3a")
    d_vb    <- tcrdist_matrix(tcrs, "human", components = "vb")
    d_combo <- tcrdist_matrix(tcrs, "human", components = c("cdr3a", "vb"))
    expect_equal(d_combo, d_cdr3a + d_vb, tolerance = 1e-10)
})

test_that("CDR3-only distances <= full distances (tcrdist_matrix)", {
    tcrs   <- make_test_tcrs()
    d_all  <- tcrdist_matrix(tcrs, "human")
    d_cdr3 <- tcrdist_matrix(tcrs, "human", components = "cdr3")
    expect_true(all(d_cdr3 <= d_all))
})


# -- tcrdist_sparse additivity tests ----------------------------------------

test_that("v_region + cdr3 = all (tcrdist_sparse, threshold=Inf)", {
    tcrs   <- make_test_tcrs()
    sp_all  <- as.matrix(tcrdist_sparse(tcrs, "human", threshold = Inf))
    sp_v    <- as.matrix(tcrdist_sparse(tcrs, "human", threshold = Inf,
                                         components = "v_region"))
    sp_cdr3 <- as.matrix(tcrdist_sparse(tcrs, "human", threshold = Inf,
                                         components = "cdr3"))
    expect_equal(sp_v + sp_cdr3, sp_all, tolerance = 1e-10)
})

test_that("alpha + beta = all (tcrdist_sparse, threshold=Inf)", {
    tcrs   <- make_test_tcrs()
    sp_all <- as.matrix(tcrdist_sparse(tcrs, "human", threshold = Inf))
    sp_a   <- as.matrix(tcrdist_sparse(tcrs, "human", threshold = Inf,
                                        components = "alpha"))
    sp_b   <- as.matrix(tcrdist_sparse(tcrs, "human", threshold = Inf,
                                        components = "beta"))
    expect_equal(sp_a + sp_b, sp_all, tolerance = 1e-10)
})


# -- tcrdist_rect additivity tests ------------------------------------------

test_that("v_region + cdr3 = all (tcrdist_rect)", {
    tcrs <- make_test_tcrs()
    ref  <- tcrs[c(1, 3), ]
    r_all  <- tcrdist_rect(tcrs, ref, "human")
    r_v    <- tcrdist_rect(tcrs, ref, "human", components = "v_region")
    r_cdr3 <- tcrdist_rect(tcrs, ref, "human", components = "cdr3")
    expect_equal(r_v + r_cdr3, r_all, tolerance = 1e-10)
})

test_that("alpha + beta = all (tcrdist_rect)", {
    tcrs <- make_test_tcrs()
    ref  <- tcrs[c(1, 3), ]
    r_all <- tcrdist_rect(tcrs, ref, "human")
    r_a   <- tcrdist_rect(tcrs, ref, "human", components = "alpha")
    r_b   <- tcrdist_rect(tcrs, ref, "human", components = "beta")
    expect_equal(r_a + r_b, r_all, tolerance = 1e-10)
})


# -- tcrdist_knn component tests --------------------------------------------

test_that("tcrdist_knn with cdr3 returns smaller distances than all", {
    tcrs     <- make_test_tcrs()
    knn_all  <- tcrdist_knn(tcrs, "human", K = 1L)
    knn_cdr3 <- tcrdist_knn(tcrs, "human", K = 1L, components = "cdr3")
    # CDR3-only KNN nearest distances should be <= full nearest distances
    expect_true(all(knn_cdr3$knn_distances <= knn_all$knn_distances))
})

test_that("tcrdist_knn components produces consistent results with tcrdist_matrix", {
    tcrs <- make_test_tcrs()
    # Full matrix for CDR3-only
    d_cdr3 <- tcrdist_matrix(tcrs, "human", components = "cdr3")
    knn_cdr3 <- tcrdist_knn(tcrs, "human", K = 2L, components = "cdr3")

    # For each row, check that KNN distances match the matrix
    for (i in seq_len(nrow(tcrs))) {
        for (k in seq_len(2L)) {
            j <- knn_cdr3$knn_indices[i, k]
            expect_equal(knn_cdr3$knn_distances[i, k], d_cdr3[i, j],
                         tolerance = 1e-10)
        }
    }
})


# -- tcrdist_radius_neighbors component tests --------------------------------

test_that("tcrdist_radius_neighbors with cdr3 finds more neighbors than all", {
    tcrs      <- make_test_tcrs()
    nbrs_all  <- tcrdist_radius_neighbors(tcrs, "human", radius = 50)
    nbrs_cdr3 <- tcrdist_radius_neighbors(tcrs, "human", radius = 50,
                                            components = "cdr3")
    # CDR3-only distances are smaller, so within the same radius there should
    # be at least as many neighbors
    for (i in seq_along(nbrs_all)) {
        expect_true(length(nbrs_cdr3[[i]]$indices) >=
                    length(nbrs_all[[i]]$indices))
    }
})

test_that("tcrdist_radius_neighbors components consistent with tcrdist_matrix", {
    tcrs   <- make_test_tcrs()
    d_cdr3 <- tcrdist_matrix(tcrs, "human", components = "cdr3")
    nbrs   <- tcrdist_radius_neighbors(tcrs, "human", radius = 200,
                                        components = "cdr3")

    for (i in seq_len(nrow(tcrs))) {
        idx <- nbrs[[i]]$indices
        dst <- nbrs[[i]]$distances
        for (k in seq_along(idx)) {
            expect_equal(dst[k], d_cdr3[i, idx[k]], tolerance = 1e-10)
        }
    }
})


# -- Cross-wrapper consistency (DASH mouse data) ----------------------------

test_that("component additivity on DASH data", {
    skip_if_not(exists("dash", where = asNamespace("tcrdistR")),
                message = "DASH dataset not available")
    data(dash, package = "tcrdistR")
    sub <- dash[1:20, ]

    d_all  <- tcrdist_matrix(sub, "mouse")
    d_v    <- tcrdist_matrix(sub, "mouse", components = "v_region")
    d_cdr3 <- tcrdist_matrix(sub, "mouse", components = "cdr3")
    expect_equal(d_v + d_cdr3, d_all, tolerance = 1e-10)

    d_a <- tcrdist_matrix(sub, "mouse", components = "alpha")
    d_b <- tcrdist_matrix(sub, "mouse", components = "beta")
    expect_equal(d_a + d_b, d_all, tolerance = 1e-10)

    d_va    <- tcrdist_matrix(sub, "mouse", components = "va")
    d_cdr3a <- tcrdist_matrix(sub, "mouse", components = "cdr3a")
    d_vb    <- tcrdist_matrix(sub, "mouse", components = "vb")
    d_cdr3b <- tcrdist_matrix(sub, "mouse", components = "cdr3b")
    expect_equal(d_va + d_cdr3a + d_vb + d_cdr3b, d_all, tolerance = 1e-10)
})
