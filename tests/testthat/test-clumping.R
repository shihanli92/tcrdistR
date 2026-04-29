# test-clumping.R — Tests for the TCR clumping pipeline (Phase 4).
#
# Tests cover:
#   - rcpp_count_nuc_matches (nucleotide prefix matching)
#   - setup_tcr_groups (chain group assignment)
#   - .single_linkage_clumping (cluster merging)
#   - rcpp_poisson_test_loop (Poisson significance testing)
#   - .get_v_cdr3_nucseq / .get_j_cdr3_nucseq (gene database lookup)
#   - .analyze_junction (junction analysis)
#   - rcpp_calc_background_distributions (background distribution computation)
#   - find_clumping (end-to-end pipeline)


# ---- Test helper: simple AA -> codon mapping --------------------------------
# Deterministic codon table for constructing nucleotide sequences from protein
# sequences in tests.  Not biologically realistic, but sufficient for
# exercising the pipeline (all positions will be N-insertions in junction
# analysis, producing breakpoints at every position for resampling).

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
# 1. rcpp_count_nuc_matches
# =============================================================================

test_that("rcpp_count_nuc_matches: exact match returns full length", {
    expect_equal(rcpp_count_nuc_matches("atcgatcg", "atcgatcg"), 8L)
})

test_that("rcpp_count_nuc_matches: mismatch at position 4 returns 4", {
    # Positions 0-3 match (a,t,c,g), position 4 mismatches (a vs t)
    expect_equal(rcpp_count_nuc_matches("atcgatcg", "atcgttcg"), 4L)
})

test_that("rcpp_count_nuc_matches: empty strings return 0", {
    expect_equal(rcpp_count_nuc_matches("", "atcg"), 0L)
    expect_equal(rcpp_count_nuc_matches("atcg", ""), 0L)
    expect_equal(rcpp_count_nuc_matches("", ""), 0L)
})

test_that("rcpp_count_nuc_matches: all mismatches returns 0", {
    expect_equal(rcpp_count_nuc_matches("aaaa", "tttt"), 0L)
})

test_that("rcpp_count_nuc_matches: custom mismatch_score changes result", {
    # "atcg" vs "axcg": match, mismatch, match, match
    # With mismatch_score = -1: scores 1, 0, 1, 2 -> best=2 at pos 4
    expect_equal(rcpp_count_nuc_matches("atcg", "axcg", -1L), 4L)
    # With mismatch_score = -4: scores 1, -3, -2, -1 -> best=1 at pos 1
    expect_equal(rcpp_count_nuc_matches("atcg", "axcg", -4L), 1L)
})

test_that("rcpp_count_nuc_matches: different lengths uses shorter", {
    expect_equal(rcpp_count_nuc_matches("atcg", "atcgatcg"), 4L)
    expect_equal(rcpp_count_nuc_matches("atcgatcg", "atcg"), 4L)
})


# =============================================================================
# 2. setup_tcr_groups
# =============================================================================

test_that("setup_tcr_groups: identical alpha chains get same agroup", {
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
    # Third has different alpha chain
    expect_false(groups$agroups[3] == groups$agroups[1])
    # All beta chains differ
    expect_equal(length(unique(groups$bgroups)), 3L)
})

test_that("setup_tcr_groups: returns 0-based indices", {
    tcr_df <- data.frame(
        va = "TRAV1-1*01", ja = "TRAJ33*01", cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV19*01",  jb = "TRBJ2-7*01", cdr3b = "CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    groups <- setup_tcr_groups(tcr_df)
    expect_equal(groups$agroups, 0L)
    expect_equal(groups$bgroups, 0L)
})

test_that("setup_tcr_groups: empty data.frame returns empty vectors", {
    tcr_df <- data.frame(
        va = character(0), ja = character(0), cdr3a = character(0),
        vb = character(0), jb = character(0), cdr3b = character(0),
        stringsAsFactors = FALSE
    )
    groups <- setup_tcr_groups(tcr_df)
    expect_equal(groups$agroups, integer(0))
    expect_equal(groups$bgroups, integer(0))
})

test_that("setup_tcr_groups: errors on missing columns", {
    tcr_df <- data.frame(va = "X", ja = "Y", stringsAsFactors = FALSE)
    expect_error(setup_tcr_groups(tcr_df), "missing required columns")
})

test_that("setup_tcr_groups: includes nucseq in grouping when present", {
    # Same protein but different nucseq -> different groups
    tcr_df <- data.frame(
        va           = c("TRAV1-1*01", "TRAV1-1*01"),
        ja           = c("TRAJ33*01",  "TRAJ33*01"),
        cdr3a        = c("CAVRDSSYKLIF", "CAVRDSSYKLIF"),
        cdr3a_nucseq = c("tgcgccatcaaa", "tgcgccatcggg"),
        vb           = c("TRBV19*01",  "TRBV19*01"),
        jb           = c("TRBJ2-7*01", "TRBJ2-7*01"),
        cdr3b        = c("CASSIRSSYEQYF", "CASSIRSSYEQYF"),
        stringsAsFactors = FALSE
    )
    groups <- setup_tcr_groups(tcr_df)
    # Different nucseq -> different agroups despite same protein
    expect_false(groups$agroups[1] == groups$agroups[2])
})


# =============================================================================
# 3. .single_linkage_clumping
# =============================================================================

test_that(".single_linkage_clumping: connected chain merges to one cluster", {
    # 5 clones: 0-1 connected, 1-2 connected (chain: 0-1-2), clone 3 isolated
    all_clumped_nbrs <- list(
        "0" = c(0L, 1L),
        "1" = c(0L, 1L, 2L),
        "2" = c(1L, 2L),
        "3" = c(3L)
    )
    is_clumped <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
    clusters <- .single_linkage_clumping(all_clumped_nbrs, 5L, is_clumped)

    # Clones 0, 1, 2 are in the same cluster
    expect_equal(clusters[1], clusters[2])
    expect_equal(clusters[2], clusters[3])
    # Clone 3 is in a different (smaller) cluster
    expect_true(clusters[4] > 0L)
    expect_false(clusters[4] == clusters[1])
    # Clone 4 is not clumped
    expect_equal(clusters[5], 0L)
    # Largest cluster (3 members) gets ID 1, smaller (1 member) gets ID 2
    expect_equal(clusters[1], 1L)
    expect_equal(clusters[4], 2L)
})

test_that(".single_linkage_clumping: isolated clumped nodes get separate clusters", {
    all_clumped_nbrs <- list(
        "0" = c(0L),
        "2" = c(2L)
    )
    is_clumped <- c(TRUE, FALSE, TRUE)
    clusters <- .single_linkage_clumping(all_clumped_nbrs, 3L, is_clumped)

    expect_true(clusters[1] > 0L)
    expect_equal(clusters[2], 0L)
    expect_true(clusters[3] > 0L)
    expect_false(clusters[1] == clusters[3])
})


# =============================================================================
# 4. rcpp_poisson_test_loop
# =============================================================================

test_that("rcpp_poisson_test_loop: detects enriched neighborhood", {
    num_clones <- 5L
    agroups <- 0:4
    bgroups <- 0:4

    # Clone 0 has 3 neighbors within distance 20
    all_nbr_indices <- list(
        as.integer(c(1, 2, 3)),
        as.integer(c(0)),
        as.integer(c(0)),
        as.integer(c(0)),
        integer(0)
    )
    all_nbr_distances <- list(
        c(10.0, 15.0, 20.0),
        c(10.0),
        c(15.0),
        c(20.0),
        numeric(0)
    )

    radii <- 24L
    # Very low background frequency -> enrichment is significant
    bg_freqs <- matrix(1e-6, nrow = num_clones, ncol = 25L)

    result <- rcpp_poisson_test_loop(
        all_nbr_indices       = all_nbr_indices,
        all_nbr_distances     = all_nbr_distances,
        bg_freqs              = bg_freqs,
        agroups               = agroups,
        bgroups               = bgroups,
        radii                 = radii,
        n_bg_pairs            = 1e10,
        pvalue_threshold      = 1.0,
        num_clones            = num_clones,
        clusters_gex_nullable = NULL,
        use_conservative_pvalues = TRUE
    )

    expect_true(is.list(result))
    expect_true(all(c("clone_index", "is_clumped", "all_raw_pvalues",
                       "pvalue_adj", "num_nbrs", "expected_num_nbrs",
                       "clump_type") %in% names(result)))

    # Clone 0 should be detected as clumped (3 nbrs with tiny bg freq)
    expect_true(result$is_clumped[1])

    # all_raw_pvalues: num_clones x num_radii
    expect_equal(dim(result$all_raw_pvalues), c(5L, 1L))

    # Clone 0's raw p-value should be very small
    expect_true(result$all_raw_pvalues[1, 1] < 0.01)
})

test_that("rcpp_poisson_test_loop: Bonferroni correction applied correctly", {
    num_clones <- 3L
    agroups <- 0:2
    bgroups <- 0:2

    all_nbr_indices <- list(
        as.integer(c(1, 2)),
        as.integer(c(0, 2)),
        as.integer(c(0, 1))
    )
    all_nbr_distances <- list(
        c(10.0, 15.0),
        c(10.0, 20.0),
        c(15.0, 20.0)
    )

    radii <- c(24L, 48L)
    bg_freqs <- matrix(1e-4, nrow = num_clones, ncol = 49L)

    result <- rcpp_poisson_test_loop(
        all_nbr_indices       = all_nbr_indices,
        all_nbr_distances     = all_nbr_distances,
        bg_freqs              = bg_freqs,
        agroups               = agroups,
        bgroups               = bgroups,
        radii                 = radii,
        n_bg_pairs            = 1e8,
        pvalue_threshold      = 1.0,
        num_clones            = num_clones,
        clusters_gex_nullable = NULL,
        use_conservative_pvalues = TRUE
    )

    # Verify Bonferroni: pval_adj = raw_pval * num_radii * num_clones
    if (length(result$clone_index) > 0L) {
        ci  <- result$clone_index[1] + 1L      # 0-based -> 1-based
        ri  <- which(radii == result$nbr_radius[1])
        raw <- result$all_raw_pvalues[ci, ri]
        expect_equal(result$pvalue_adj[1],
                     raw * length(radii) * num_clones,
                     tolerance = 1e-10)
    }
})

test_that("rcpp_poisson_test_loop: no clumps when background is high", {
    num_clones <- 3L
    agroups <- 0:2
    bgroups <- 0:2

    all_nbr_indices   <- list(as.integer(c(1)), as.integer(c(0)), integer(0))
    all_nbr_distances <- list(c(10.0), c(10.0), numeric(0))

    radii <- 24L
    # High background -> 1 neighbor is expected, not enriched
    bg_freqs <- matrix(0.9, nrow = num_clones, ncol = 25L)

    result <- rcpp_poisson_test_loop(
        all_nbr_indices       = all_nbr_indices,
        all_nbr_distances     = all_nbr_distances,
        bg_freqs              = bg_freqs,
        agroups               = agroups,
        bgroups               = bgroups,
        radii                 = radii,
        n_bg_pairs            = 1e6,
        pvalue_threshold      = 0.05,
        num_clones            = num_clones,
        clusters_gex_nullable = NULL,
        use_conservative_pvalues = TRUE
    )

    # High background: no clumps should be significant at 0.05
    expect_equal(length(result$clone_index), 0L)
    expect_true(all(!result$is_clumped))
})


# =============================================================================
# 5. Gene database integration tests
# =============================================================================

test_that(".get_v_cdr3_nucseq: known human V gene returns non-empty lowercase", {
    nucseq <- .get_v_cdr3_nucseq("human", "TRAV1-1*01")
    expect_true(nchar(nucseq) > 0L)
    expect_true(grepl("^[acgt]+$", nucseq))
})

test_that(".get_j_cdr3_nucseq: known human J gene returns non-empty lowercase", {
    nucseq <- .get_j_cdr3_nucseq("human", "TRAJ33*01")
    expect_true(nchar(nucseq) > 0L)
    expect_true(grepl("^[acgt]+$", nucseq))
})

test_that(".get_v_cdr3_nucseq: unknown gene returns empty string", {
    expect_equal(.get_v_cdr3_nucseq("human", "FAKE_GENE*01"), "")
})

test_that(".get_j_cdr3_nucseq: unknown gene returns empty string", {
    expect_equal(.get_j_cdr3_nucseq("human", "FAKE_GENE*01"), "")
})


# =============================================================================
# 6. Junction analysis integration tests
# =============================================================================

test_that(".analyze_junction: alpha chain returns correct structure", {
    nucseq <- .test_prot_to_nuc("CAVRDSSYKLIF")
    result <- .analyze_junction("human", "TRAV1-1*01", "TRAJ33*01",
                                 "CAVRDSSYKLIF", nucseq)

    expect_true(is.list(result))
    expect_true(all(c("trims", "inserts", "cdr3_nucseq_src",
                       "new_nucseq", "cdr3_protseq_masked") %in% names(result)))
    expect_equal(length(result$trims), 4L)
    expect_equal(length(result$inserts), 4L)
    # Source annotation length == nucseq length
    expect_equal(nchar(result$cdr3_nucseq_src), nchar(nucseq))
    # Source chars are V, J, N, or D
    expect_true(grepl("^[VJND]+$", result$cdr3_nucseq_src))
})

test_that(".analyze_junction: beta chain returns correct structure", {
    nucseq <- .test_prot_to_nuc("CASSIRSSYEQYF")
    result <- .analyze_junction("human", "TRBV19*01", "TRBJ2-7*01",
                                 "CASSIRSSYEQYF", nucseq)

    expect_true(is.list(result))
    expect_equal(length(result$trims), 4L)
    expect_equal(length(result$inserts), 4L)
    expect_equal(nchar(result$cdr3_nucseq_src), nchar(nucseq))
    expect_true(grepl("^[VJND]+$", result$cdr3_nucseq_src))
})

test_that(".analyze_junction: V/J trims are non-negative", {
    nucseq <- .test_prot_to_nuc("CAVRDSSYKLIF")
    result <- .analyze_junction("human", "TRAV1-1*01", "TRAJ33*01",
                                 "CAVRDSSYKLIF", nucseq)
    # trims: c(v_trim, d0_trim, d1_trim, j_trim) -- all non-negative
    expect_true(all(result$trims >= 0L))
})


# =============================================================================
# 7. Background distribution computation
# =============================================================================

test_that("rcpp_calc_background_distributions: correct shape and monotonicity", {
    v_dist_a <- .compute_v_region_distance_matrix("human", "A")
    v_dist_b <- .compute_v_region_distance_matrix("human", "B")

    fg_va    <- c("TRAV1-1*01", "TRAV12-2*01")
    fg_cdr3a <- c("CAVRDSSYKLIF", "CAVSANSGTYF")
    fg_vb    <- c("TRBV19*01", "TRBV20-1*01")
    fg_cdr3b <- c("CASSIRSSYEQYF", "CSARDRTGNTIYF")

    # Small background: repeat 2 chains 10 times each (20 total per chain type)
    bg_va    <- rep(c("TRAV1-1*01", "TRAV12-2*01"), 10L)
    bg_cdr3a <- rep(c("CAVRDSSYKLIF", "CAVSANSGTYF"), 10L)
    bg_vb    <- rep(c("TRBV19*01", "TRBV20-1*01"), 10L)
    bg_cdr3b <- rep(c("CASSIRSSYEQYF", "CSARDRTGNTIYF"), 10L)

    max_dist <- 48L

    result <- rcpp_calc_background_distributions(
        fg_va, fg_cdr3a, fg_vb, fg_cdr3b,
        bg_va, bg_cdr3a, bg_vb, bg_cdr3b,
        v_dist_a, v_dist_b,
        max_dist
    )

    # Shape: n_fg x (max_dist + 1)
    expect_equal(dim(result), c(2L, 49L))

    # Values should be non-decreasing across columns (cumulative frequencies)
    for (i in seq_len(nrow(result))) {
        diffs <- diff(result[i, ])
        expect_true(all(diffs >= -1e-10),
                    info = sprintf("Row %d not monotonically non-decreasing", i))
    }

    # All values in [0, 1]
    expect_true(all(result >= 0))
    expect_true(all(result <= 1 + 1e-10))

    # First column (dist=0) should be positive (from pseudocount) but small
    expect_true(all(result[, 1] > 0))
    expect_true(all(result[, 1] < 0.5))
})


# =============================================================================
# 8. find_clumping end-to-end
# =============================================================================

test_that("find_clumping: runs to completion and returns correct structure", {
    skip_on_cran()
    set.seed(42L)

    # ---- Build 8 TCRs: 5 tight cluster + 3 diverse ----
    # Cluster: same V/J genes, CDR3s differ at 1 position each.
    # Each clone has a unique alpha AND beta CDR3 so that
    # setup_tcr_groups assigns different groups (no same-group masking).
    cluster_n  <- 5L
    diverse_n  <- 3L
    total_n    <- cluster_n + diverse_n

    cluster_va    <- rep("TRAV1-1*01",  cluster_n)
    cluster_ja    <- rep("TRAJ33*01",   cluster_n)
    cluster_vb    <- rep("TRBV19*01",   cluster_n)
    cluster_jb    <- rep("TRBJ2-7*01",  cluster_n)
    cluster_cdr3a <- c("CAVRDSSYKLIF", "CAVRDTSYKLIF", "CAVRDASYKLIF",
                        "CAVRDGSYKLIF", "CAVRDNSYKLIF")
    cluster_cdr3b <- c("CASSIRSSYEQYF", "CASSIRTSYEQYF", "CASSIRASYEQYF",
                        "CASSIRGSYEQYF", "CASSIRNSYEQYF")

    diverse_va    <- rep("TRAV12-2*01", diverse_n)
    diverse_ja    <- rep("TRAJ49*01",   diverse_n)
    diverse_vb    <- rep("TRBV20-1*01", diverse_n)
    diverse_jb    <- rep("TRBJ1-1*01",  diverse_n)
    diverse_cdr3a <- c("CAVSANSGTYF",  "CAVSNFRGTYF",  "CAVSQDRGTYF")
    diverse_cdr3b <- c("CSARDRTGNTIYF", "CSARDHTNNTIYF", "CSARDPFGNTIYF")

    all_cdr3a <- c(cluster_cdr3a, diverse_cdr3a)
    all_cdr3b <- c(cluster_cdr3b, diverse_cdr3b)

    tcr_df <- data.frame(
        va           = c(cluster_va, diverse_va),
        ja           = c(cluster_ja, diverse_ja),
        cdr3a        = all_cdr3a,
        cdr3a_nucseq = vapply(all_cdr3a, .test_prot_to_nuc, character(1L)),
        vb           = c(cluster_vb, diverse_vb),
        jb           = c(cluster_jb, diverse_jb),
        cdr3b        = all_cdr3b,
        cdr3b_nucseq = vapply(all_cdr3b, .test_prot_to_nuc, character(1L)),
        stringsAsFactors = FALSE
    )
    rownames(tcr_df) <- NULL

    radii <- c(48L, 96L)

    result <- find_clumping(
        tcr_df, "human",
        radii              = radii,
        num_random_samples = 500L,
        pvalue_threshold   = 1.0,
        verbose            = FALSE
    )

    # ---- Structure checks ----
    expect_true(is.list(result))
    expect_true(all(c("results_df", "is_clumped", "clusters",
                       "all_raw_pvalues") %in% names(result)))

    expect_equal(length(result$is_clumped), total_n)
    expect_equal(length(result$clusters), total_n)
    expect_equal(nrow(result$all_raw_pvalues), total_n)
    expect_equal(ncol(result$all_raw_pvalues), length(radii))

    expected_cols <- c("clump_type", "clone_index", "nbr_radius", "pvalue_adj",
                       "num_nbrs", "expected_num_nbrs", "raw_count",
                       "va", "ja", "cdr3a", "vb", "jb", "cdr3b",
                       "clumping_group", "clonotype_fdr_value")
    expect_true(all(expected_cols %in% colnames(result$results_df)))

    # Clusters are non-negative integers
    expect_true(all(result$clusters >= 0L))
    # is_clumped is logical of correct length
    expect_true(is.logical(result$is_clumped))
})

test_that("find_clumping: errors on missing nucleotide columns", {
    tcr_df <- data.frame(
        va = "TRAV1-1*01", ja = "TRAJ33*01", cdr3a = "CAVRDSSYKLIF",
        vb = "TRBV19*01",  jb = "TRBJ2-7*01", cdr3b = "CASSIRSSYEQYF",
        stringsAsFactors = FALSE
    )
    expect_error(find_clumping(tcr_df, "human"), "missing required columns")
})

test_that("find_clumping: empty data.frame returns empty results", {
    tcr_df <- data.frame(
        va = character(0), ja = character(0), cdr3a = character(0),
        cdr3a_nucseq = character(0),
        vb = character(0), jb = character(0), cdr3b = character(0),
        cdr3b_nucseq = character(0),
        stringsAsFactors = FALSE
    )
    result <- find_clumping(tcr_df, "human", verbose = FALSE)
    expect_equal(nrow(result$results_df), 0L)
    expect_equal(length(result$is_clumped), 0L)
    expect_equal(length(result$clusters), 0L)
})
