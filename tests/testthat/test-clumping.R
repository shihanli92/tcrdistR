# test-clumping.R — Tests for the TCR clumping pipeline.
# Internal helpers (.analyze_junction, .get_v_cdr3_nucseq, .get_j_cdr3_nucseq)
# are tested implicitly through find_clumping.


# ---- Test helper: simple AA -> codon mapping --------------------------------

.test_aa_to_codon <- c(
    A = "gcc", C = "tgc", D = "gac", E = "gag", F = "ttc", G = "ggc",
    H = "cac", I = "atc", K = "aag", L = "ctg", M = "atg", N = "aac",
    P = "ccc", Q = "cag", R = "cgg", S = "agc", T = "acc", V = "gtc",
    W = "tgg", Y = "tac"
)

.test_prot_to_nuc <- function(protseq) {
    aas <- strsplit(protseq, "", fixed = TRUE)[[1L]]
    paste0(.test_aa_to_codon[aas], collapse = "")
}


# =============================================================================
# rcpp_count_nuc_matches
# =============================================================================

test_that("rcpp_count_nuc_matches handles all match scenarios", {
    expect_equal(rcpp_count_nuc_matches("atcgatcg", "atcgatcg"), 8L)
    expect_equal(rcpp_count_nuc_matches("atcgatcg", "atcgttcg"), 4L)
    expect_equal(rcpp_count_nuc_matches("", "atcg"), 0L)
    expect_equal(rcpp_count_nuc_matches("atcg", ""), 0L)
    expect_equal(rcpp_count_nuc_matches("", ""), 0L)
    expect_equal(rcpp_count_nuc_matches("aaaa", "tttt"), 0L)
    # Different lengths: uses shorter
    expect_equal(rcpp_count_nuc_matches("atcg", "atcgatcg"), 4L)
    # Custom mismatch_score
    expect_equal(rcpp_count_nuc_matches("atcg", "axcg", -1L), 4L)
    expect_equal(rcpp_count_nuc_matches("atcg", "axcg", -4L), 1L)
})


# =============================================================================
# setup_tcr_groups
# =============================================================================

test_that("setup_tcr_groups assigns groups correctly", {
    tcr_df <- data.frame(
        va    = c("TRAV1-1*01", "TRAV1-1*01", "TRAV12-2*01"),
        ja    = c("TRAJ33*01",  "TRAJ33*01",  "TRAJ49*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF", "CAVSANSGTYF"),
        vb    = c("TRBV19*01",   "TRBV20-1*01",  "TRBV19*01"),
        jb    = c("TRBJ2-7*01",  "TRBJ1-1*01",   "TRBJ2-7*01"),
        cdr3b = c("CASSIRSSYEQYF", "CSARDRTGNTIYF", "CASSLGQAYEQYF"),
        stringsAsFactors = FALSE
    )
    groups <- setup_tcr_groups(tcr_df)

    # First two share identical alpha chains
    expect_equal(groups$agroups[1], groups$agroups[2])
    expect_false(groups$agroups[3] == groups$agroups[1])
    expect_equal(length(unique(groups$bgroups)), 3L)

    # 0-based indices
    single <- data.frame(
        va = "TRAV1-1*01", ja = "TRAJ33*01", cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV19*01",  jb = "TRBJ2-7*01", cdr3b = "CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    expect_equal(setup_tcr_groups(single)$agroups, 0L)

    # Nucseq affects grouping
    tcr_nuc <- data.frame(
        va = c("TRAV1-1*01", "TRAV1-1*01"),
        ja = c("TRAJ33*01", "TRAJ33*01"),
        cdr3a = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        cdr3a_nucseq = c("tgcgccatcaaa", "tgcgccatcggg"),
        vb = c("TRBV19*01", "TRBV19*01"),
        jb = c("TRBJ2-7*01", "TRBJ2-7*01"),
        cdr3b = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )
    expect_false(setup_tcr_groups(tcr_nuc)$agroups[1] ==
                     setup_tcr_groups(tcr_nuc)$agroups[2])

    # Missing columns error
    expect_error(setup_tcr_groups(data.frame(va = "X", ja = "Y",
                                              stringsAsFactors = FALSE)),
                 "missing required columns")
})


# =============================================================================
# .single_linkage_clumping
# =============================================================================

test_that(".single_linkage_clumping merges connected components correctly", {
    # Chain: 0-1-2 connected, clone 3 isolated, clone 4 not clumped
    all_clumped_nbrs <- list(
        "0" = c(0L, 1L), "1" = c(0L, 1L, 2L),
        "2" = c(1L, 2L), "3" = c(3L)
    )
    is_clumped <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
    clusters <- .single_linkage_clumping(all_clumped_nbrs, 5L, is_clumped)

    expect_equal(clusters[1], clusters[2])
    expect_equal(clusters[2], clusters[3])
    expect_true(clusters[4] > 0L)
    expect_false(clusters[4] == clusters[1])
    expect_equal(clusters[5], 0L)
    # Largest cluster gets ID 1
    expect_equal(clusters[1], 1L)
    expect_equal(clusters[4], 2L)

    # Isolated clumped nodes
    nbrs2 <- list("0" = c(0L), "2" = c(2L))
    cl2 <- .single_linkage_clumping(nbrs2, 3L, c(TRUE, FALSE, TRUE))
    expect_true(cl2[1] > 0L && cl2[3] > 0L)
    expect_equal(cl2[2], 0L)
    expect_false(cl2[1] == cl2[3])
})


# =============================================================================
# rcpp_poisson_test_loop
# =============================================================================

test_that("rcpp_poisson_test_loop detects enrichment and applies Bonferroni", {
    num_clones <- 5L
    agroups <- 0:4
    bgroups <- 0:4

    all_nbr_indices <- list(
        as.integer(c(1, 2, 3)), as.integer(c(0)),
        as.integer(c(0)), as.integer(c(0)), integer(0)
    )
    all_nbr_distances <- list(
        c(10.0, 15.0, 20.0), c(10.0), c(15.0), c(20.0), numeric(0)
    )

    bg_freqs <- matrix(1e-6, nrow = num_clones, ncol = 25L)

    result <- rcpp_poisson_test_loop(
        all_nbr_indices = all_nbr_indices,
        all_nbr_distances = all_nbr_distances,
        bg_freqs = bg_freqs, agroups = agroups, bgroups = bgroups,
        radii = 24L, n_bg_pairs = 1e10, pvalue_threshold = 1.0,
        num_clones = num_clones, clusters_gex_nullable = NULL,
        use_conservative_pvalues = TRUE
    )

    expect_true(is.list(result))
    expect_true(result$is_clumped[1])
    expect_equal(dim(result$all_raw_pvalues), c(5L, 1L))
    expect_true(result$all_raw_pvalues[1, 1] < 0.01)

    # Multi-radius Bonferroni
    radii2 <- c(24L, 48L)
    bg2 <- matrix(1e-4, nrow = 3L, ncol = 49L)
    nbr_idx2 <- list(as.integer(c(1, 2)), as.integer(c(0, 2)),
                      as.integer(c(0, 1)))
    nbr_dst2 <- list(c(10.0, 15.0), c(10.0, 20.0), c(15.0, 20.0))

    r2 <- rcpp_poisson_test_loop(
        all_nbr_indices = nbr_idx2, all_nbr_distances = nbr_dst2,
        bg_freqs = bg2, agroups = 0:2, bgroups = 0:2,
        radii = radii2, n_bg_pairs = 1e8, pvalue_threshold = 1.0,
        num_clones = 3L, clusters_gex_nullable = NULL,
        use_conservative_pvalues = TRUE
    )
    if (length(r2$clone_index) > 0L) {
        ci <- r2$clone_index[1] + 1L
        ri <- which(radii2 == r2$nbr_radius[1])
        raw <- r2$all_raw_pvalues[ci, ri]
        expect_equal(r2$pvalue_adj[1], raw * length(radii2) * 3L,
                     tolerance = 1e-10)
    }
})

test_that("rcpp_poisson_test_loop: no clumps when background is high", {
    nbr_idx <- list(as.integer(c(1)), as.integer(c(0)), integer(0))
    nbr_dst <- list(c(10.0), c(10.0), numeric(0))
    bg <- matrix(0.9, nrow = 3L, ncol = 25L)

    result <- rcpp_poisson_test_loop(
        all_nbr_indices = nbr_idx, all_nbr_distances = nbr_dst,
        bg_freqs = bg, agroups = 0:2, bgroups = 0:2,
        radii = 24L, n_bg_pairs = 1e6, pvalue_threshold = 0.05,
        num_clones = 3L, clusters_gex_nullable = NULL,
        use_conservative_pvalues = TRUE
    )
    expect_equal(length(result$clone_index), 0L)
    expect_true(all(!result$is_clumped))
})


# =============================================================================
# rcpp_calc_background_distributions
# =============================================================================

test_that("rcpp_calc_background_distributions: correct shape and monotonicity", {
    v_dist_a <- .compute_v_region_distance_matrix("human", "A")
    v_dist_b <- .compute_v_region_distance_matrix("human", "B")

    fg_va    <- c("TRAV1-1*01", "TRAV12-2*01")
    fg_cdr3a <- c("CAVRDSSYKLIF", "CAVSANSGTYF")
    fg_vb    <- c("TRBV19*01", "TRBV20-1*01")
    fg_cdr3b <- c("CASSIRSSYEQYF", "CSARDRTGNTIYF")

    bg_va    <- rep(c("TRAV1-1*01", "TRAV12-2*01"), 10L)
    bg_cdr3a <- rep(c("CAVRDSSYKLIF", "CAVSANSGTYF"), 10L)
    bg_vb    <- rep(c("TRBV19*01", "TRBV20-1*01"), 10L)
    bg_cdr3b <- rep(c("CASSIRSSYEQYF", "CSARDRTGNTIYF"), 10L)

    result <- rcpp_calc_background_distributions(
        fg_va, fg_cdr3a, fg_vb, fg_cdr3b,
        bg_va, bg_cdr3a, bg_vb, bg_cdr3b,
        v_dist_a, v_dist_b, 48L
    )

    expect_equal(dim(result), c(2L, 49L))
    for (i in seq_len(nrow(result))) {
        expect_true(all(diff(result[i, ]) >= -1e-10))
    }
    expect_true(all(result >= 0 & result <= 1 + 1e-10))
    expect_true(all(result[, 1] > 0 & result[, 1] < 0.5))
})


# =============================================================================
# find_clumping end-to-end
# =============================================================================

test_that("find_clumping: runs to completion and returns correct structure", {
    skip_on_cran()
    set.seed(42L)

    cluster_n <- 5L
    diverse_n <- 3L
    total_n <- cluster_n + diverse_n

    cluster_cdr3a <- c("CAVRDSSYKLIF", "CAVRDTSYKLIF", "CAVRDASYKLIF",
                        "CAVRDGSYKLIF", "CAVRDNSYKLIF")
    cluster_cdr3b <- c("CASSIRSSYEQYF", "CASSIRTSYEQYF", "CASSIRASYEQYF",
                        "CASSIRGSYEQYF", "CASSIRNSYEQYF")
    diverse_cdr3a <- c("CAVSANSGTYF", "CAVSNFRGTYF", "CAVSQDRGTYF")
    diverse_cdr3b <- c("CSARDRTGNTIYF", "CSARDHTNNTIYF", "CSARDPFGNTIYF")

    all_cdr3a <- c(cluster_cdr3a, diverse_cdr3a)
    all_cdr3b <- c(cluster_cdr3b, diverse_cdr3b)

    tcr_df <- data.frame(
        va = c(rep("TRAV1-1*01", cluster_n), rep("TRAV12-2*01", diverse_n)),
        ja = c(rep("TRAJ33*01", cluster_n), rep("TRAJ49*01", diverse_n)),
        cdr3a = all_cdr3a,
        cdr3a_nucseq = vapply(all_cdr3a, .test_prot_to_nuc, character(1L)),
        vb = c(rep("TRBV19*01", cluster_n), rep("TRBV20-1*01", diverse_n)),
        jb = c(rep("TRBJ2-7*01", cluster_n), rep("TRBJ1-1*01", diverse_n)),
        cdr3b = all_cdr3b,
        cdr3b_nucseq = vapply(all_cdr3b, .test_prot_to_nuc, character(1L)),
        stringsAsFactors = FALSE
    )
    rownames(tcr_df) <- NULL

    result <- find_clumping(
        tcr_df, "human", radii = c(48L, 96L),
        num_random_samples = 500L, pvalue_threshold = 1.0, verbose = FALSE
    )

    expect_true(is.list(result))
    expect_equal(length(result$is_clumped), total_n)
    expect_equal(length(result$clusters), total_n)
    expect_equal(nrow(result$all_raw_pvalues), total_n)
    expect_equal(ncol(result$all_raw_pvalues), 2L)
    expect_true(all(result$clusters >= 0L))
    expect_true(is.logical(result$is_clumped))

    expected_cols <- c("clump_type", "clone_index", "nbr_radius", "pvalue_adj",
                       "num_nbrs", "expected_num_nbrs", "raw_count",
                       "va", "ja", "cdr3a", "vb", "jb", "cdr3b",
                       "clumping_group", "clonotype_fdr_value")
    expect_true(all(expected_cols %in% colnames(result$results_df)))
})

test_that("find_clumping: error and edge cases", {
    # Missing nucseq columns
    tcr_df_no_nuc <- data.frame(
        va = "TRAV1-1*01", ja = "TRAJ33*01", cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV19*01", jb = "TRBJ2-7*01", cdr3b = "CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    expect_error(find_clumping(tcr_df_no_nuc, "human"),
                 "missing required columns")

    # Empty data.frame
    tcr_df_empty <- data.frame(
        va = character(0), ja = character(0), cdr3a = character(0),
        cdr3a_nucseq = character(0),
        vb = character(0), jb = character(0), cdr3b = character(0),
        cdr3b_nucseq = character(0),
        stringsAsFactors = FALSE
    )
    result_empty <- find_clumping(tcr_df_empty, "human", verbose = FALSE)
    expect_equal(nrow(result_empty$results_df), 0L)
    expect_equal(length(result_empty$is_clumped), 0L)
})
