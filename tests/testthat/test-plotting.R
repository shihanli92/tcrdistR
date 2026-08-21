# Tests for visualization functions.
# Tests verify structure and correctness, not visual appearance.
# BLOSUM62 matrix tests live in test-blosum.R (not duplicated here).


# ===========================================================================
# Internal helpers: palette, alignment, PWM
# ===========================================================================

test_that(".tcrdistR_palette returns valid colors with correct tab switching", {
    pal5  <- tcrdistR:::.tcrdistR_palette(5)
    pal10 <- tcrdistR:::.tcrdistR_palette(10)
    pal15 <- tcrdistR:::.tcrdistR_palette(15)
    expect_length(pal5, 5L)
    expect_length(pal10, 10L)
    expect_length(pal15, 15L)
    expect_true(all(grepl("^#[0-9a-fA-F]{6}$", pal15)))
    # tab10 first color
    expect_equal(pal10[1L], "#1f77b4")
    # tab20 second color (lighter variant)
    pal11 <- tcrdistR:::.tcrdistR_palette(11)
    expect_equal(pal11[2L], "#aec7e8")
})

test_that(".align_cdr3_regions handles equal, unequal, and custom gap", {
    # Equal length: unchanged
    eq <- tcrdistR:::.align_cdr3_regions("ABCDE", "FGHIJ")
    expect_equal(eq$a, "ABCDE")
    expect_equal(eq$b, "FGHIJ")
    # Unequal: gaps in shorter
    uneq <- tcrdistR:::.align_cdr3_regions("CASSI", "CASSILY")
    expect_equal(nchar(uneq$a), nchar(uneq$b))
    expect_true(grepl("-", uneq$a))
    expect_false(grepl("-", uneq$b))
    n_gaps <- nchar(gsub("[^-]", "", uneq$a))
    expect_equal(n_gaps, 2L)
    # Single char shorter
    sc <- tcrdistR:::.align_cdr3_regions("A", "ABC")
    expect_equal(nchar(sc$a), 3L)
    # Custom gap character
    cg <- tcrdistR:::.align_cdr3_regions("CASSI", "CASSILY",
                                          gap_character = ".")
    expect_true(grepl("\\.", cg$a))
})

test_that(".build_cdr3_pwm produces valid PWM with trim and alignment", {
    # Single sequence: identity PWM
    r1 <- tcrdistR:::.build_cdr3_pwm("CASSIRSSYEQY", trim = FALSE)
    expect_true(is.matrix(r1$pwm))
    expect_equal(r1$center_idx, 1L)
    for (j in seq_len(ncol(r1$pwm))) {
        expect_equal(sum(r1$pwm[, j]), 1.0, tolerance = 1e-10)
    }
    # Multiple sequences: columns sum to 1, aligned_seqs have equal length
    seqs <- c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY")
    r2 <- tcrdistR:::.build_cdr3_pwm(seqs, trim = FALSE)
    expect_true(all(abs(colSums(r2$pwm) - 1.0) < 1e-10))
    expect_true(r2$center_idx >= 1L && r2$center_idx <= length(seqs))
    expect_true(length(unique(nchar(r2$aligned_seqs))) == 1L)
    # Trim removes 3+2 residues
    seq_full <- "CASSIRSSYEQYF"
    r_trim <- tcrdistR:::.build_cdr3_pwm(seq_full, trim = TRUE)
    r_notrim <- tcrdistR:::.build_cdr3_pwm(seq_full, trim = FALSE)
    expect_equal(ncol(r_trim$pwm), nchar(seq_full) - 5L)
    expect_equal(ncol(r_notrim$pwm), nchar(seq_full))
})


# ===========================================================================
# Exported plot functions
# ===========================================================================

test_that("plot_cdr3_length handles alpha, both chains, and empty input", {
    skip_if_not_installed("ggplot2")

    tcr_df <- data.frame(
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY"),
        stringsAsFactors = FALSE
    )
    expect_s3_class(plot_cdr3_length(tcr_df, chain = "alpha"), "gg")
    expect_s3_class(plot_cdr3_length(tcr_df, chain = "both"), "gg")

    empty_df <- data.frame(cdr3a = character(0), cdr3b = character(0),
                            stringsAsFactors = FALSE)
    expect_s3_class(plot_cdr3_length(empty_df, chain = "both"), "gg")
})

test_that("plot_gene_usage handles normal, empty, and max_genes", {
    skip_if_not_installed("ggplot2")

    tcr_df <- data.frame(
        va = c("TRAV1-2*01", "TRAV1-2*02", "TRAV12-1*01",
               "TRAV1-2*01", "TRAV12-1*01"),
        stringsAsFactors = FALSE
    )
    expect_s3_class(plot_gene_usage(tcr_df, "va"), "gg")

    empty_df <- data.frame(va = character(0), stringsAsFactors = FALSE)
    expect_s3_class(plot_gene_usage(empty_df, "va"), "gg")

    genes <- paste0("TRAV", seq_len(30))
    big_df <- data.frame(va = genes, stringsAsFactors = FALSE)
    p <- plot_gene_usage(big_df, "va", strip_allele = FALSE, max_genes = 5L)
    pb <- ggplot2::ggplot_build(p)
    expect_true(nrow(pb$data[[1L]]) <= 5L)
})

test_that("plot_tcrdist_heatmap works with and without clustering", {
    skip_if_not_installed("ggplot2")

    mat <- matrix(c(0, 10, 20, 10, 0, 15, 20, 15, 0), nrow = 3)
    expect_s3_class(plot_tcrdist_heatmap(mat), "gg")
    expect_s3_class(plot_tcrdist_heatmap(mat, cluster = FALSE,
                                          labels = c("A", "B", "C")), "gg")
})

test_that("plot_distance_distribution handles params and vector input", {
    skip_if_not_installed("ggplot2")

    mat <- matrix(c(0, 10, 20, 30, 10, 0, 15, 25,
                     20, 15, 0, 12, 30, 25, 12, 0), nrow = 4)
    expect_s3_class(plot_distance_distribution(mat), "gg")
    expect_s3_class(plot_distance_distribution(mat, binwidth = 5), "gg")
    p <- plot_distance_distribution(mat, threshold = 18)
    expect_s3_class(p, "gg")
    expect_true(length(p$layers) > 3L)

    dists <- c(5, 10, 15, 20, 25, 30)
    expect_s3_class(plot_distance_distribution(dists, threshold = 12), "gg")
})

test_that("plot_tcr_scatter handles coloring modes", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    # No color, continuous, categorical, custom palette
    expect_s3_class(plot_tcr_scatter(coords), "gg")
    expect_s3_class(plot_tcr_scatter(coords, color_by = runif(10)), "gg")
    expect_s3_class(plot_tcr_scatter(coords,
                     color_by = rep(c("A", "B"), 5)), "gg")
    expect_s3_class(plot_tcr_scatter(coords,
                     color_by = rep(c("A", "B"), 5),
                     palette = c("red", "blue")), "gg")
})

test_that("plot_tcr_scatter handles faceting and highlighting", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    # Facet by vector
    expect_s3_class(
        plot_tcr_scatter(coords, facet_by = rep(c("A", "B"), 5)), "gg")
    # Facet by 2-column data.frame
    facets <- data.frame(row = rep(c("R1", "R2"), 5),
                          col = rep(c("C1", "C2"), each = 5))
    expect_s3_class(
        plot_tcr_scatter(coords, color_by = rep(c("X", "Y"), 5),
                          facet_by = facets), "gg")
    # Highlight: logical, integer, with color_by
    hl <- c(TRUE, TRUE, TRUE, rep(FALSE, 7))
    expect_s3_class(plot_tcr_scatter(coords, highlight = hl), "gg")
    expect_s3_class(plot_tcr_scatter(coords, highlight = c(1L, 5L, 10L)), "gg")
    expect_s3_class(
        plot_tcr_scatter(coords, color_by = rep(c("A", "B"), 5),
                          highlight = hl), "gg")
})

test_that("plot_tcr_scatter metadata lookup and named-list highlight", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    meta <- data.frame(group = rep(c("A", "B"), 5),
                        site = rep(c("X", "Y", "X", "Y", "X"), 2),
                        stringsAsFactors = FALSE)
    expect_s3_class(
        plot_tcr_scatter(coords, color_by = "group", facet_by = "group",
                          metadata = meta), "gg")
    expect_s3_class(
        plot_tcr_scatter(coords, highlight = list(group = "A", site = "X"),
                          metadata = meta), "gg")
})


# ===========================================================================
# CDR3 logo (requires ggseqlogo)
# ===========================================================================

test_that("plot_cdr3_logo works with bits, prob, and single sequence", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")

    seqs <- c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF",
              "CAVRDSSYKLIF", "CAVKDSYKLIF")
    expect_s3_class(plot_cdr3_logo(seqs, chain = "alpha", method = "bits"),
                    "gg")
    expect_s3_class(plot_cdr3_logo(seqs, chain = "beta", method = "prob"),
                    "gg")
    expect_s3_class(plot_cdr3_logo("CAVRDSSYKLIF", chain = "alpha"), "gg")
})


# ===========================================================================
# Dendrogram
# ===========================================================================

test_that("plot_tcrdist_dendrogram works with and without color_by", {
    skip_if_not_installed("ggplot2")
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

    expect_s3_class(plot_tcrdist_dendrogram(tcr_df, "human"), "gg")
    expect_s3_class(
        plot_tcrdist_dendrogram(tcr_df, "human",
                                 color_by = c("A", "A", "B", "B", "A")), "gg")
})


# ===========================================================================
# Integration: kernel PCA -> scatter
# ===========================================================================

test_that("kernel PCA to scatter pipeline works", {
    skip_if_not_installed("ggplot2")
    skip_on_cran()

    fixture_path <- test_path("fixtures", "dash.csv")
    skip_if(!file.exists(fixture_path), "dash.csv fixture not found")

    dash <- utils::read.csv(fixture_path, stringsAsFactors = FALSE,
                             nrows = 20)
    tcr_df <- data.frame(
        va    = dash$v_a_gene,
        cdr3a = dash$cdr3_a_aa,
        vb    = dash$v_b_gene,
        cdr3b = dash$cdr3_b_aa,
        stringsAsFactors = FALSE
    )

    kpca <- compute_tcrdist_kernel_pca(tcr_df, "mouse",
                                        n_components = 5L)
    p <- plot_tcr_scatter(kpca$embeddings[, 1:2],
                           color_by = dash$epitope,
                           title = "Kernel PCA")
    expect_s3_class(p, "gg")
})


# ===========================================================================
# plot_vj_gene_logo
# ===========================================================================

test_that("plot_vj_gene_logo works for V/J genes with max_genes and organisms", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("png")
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ]

    expect_s3_class(
        plot_vj_gene_logo(pa$vb, organism = "mouse",
                           gene_type = "V", chain = "beta"), "gg")
    expect_s3_class(
        plot_vj_gene_logo(pa$ja, organism = "mouse",
                           gene_type = "J", chain = "alpha"), "gg")
    expect_s3_class(
        plot_vj_gene_logo(dash$vb, organism = "mouse",
                           gene_type = "V", chain = "beta",
                           max_genes = 3L), "gg")
    # Human organism
    genes <- rep(c("TRAV1-2*01", "TRAV12-1*01", "TRAV10*01"), c(5, 3, 2))
    expect_s3_class(
        plot_vj_gene_logo(genes, organism = "human",
                           gene_type = "V", chain = "alpha"), "gg")
})


# ===========================================================================
# compute_nucseq_src
# ===========================================================================

test_that("compute_nucseq_src returns correct labels for alpha and beta", {
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:5, ]

    src_a <- compute_nucseq_src(pa, organism = "mouse", chain = "alpha")
    expect_type(src_a, "list")
    expect_length(src_a, 5L)
    for (s in src_a) {
        if (!is.null(s)) expect_true(all(s %in% c("V", "N", "J")))
    }

    src_b <- compute_nucseq_src(pa, organism = "mouse", chain = "beta")
    for (s in src_b) {
        if (!is.null(s)) expect_true(all(s %in% c("V", "N1", "D", "N2", "J")))
    }

    # Missing nucseq returns NULL
    pa$cdr3b_nucseq <- NA_character_
    src_na <- compute_nucseq_src(pa, organism = "mouse", chain = "beta")
    expect_true(all(vapply(src_na, is.null, logical(1L))))
})


# ===========================================================================
# plot_cdr3_logo with return_junction_pwm
# ===========================================================================

test_that("plot_cdr3_logo return_junction_pwm works and FALSE is compatible", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:10, ]
    src <- compute_nucseq_src(pa, organism = "mouse", chain = "beta")

    result <- plot_cdr3_logo(pa$cdr3b, chain = "beta",
                              nucseq_src = src,
                              return_junction_pwm = TRUE)
    expect_type(result, "list")
    expect_true(all(c("plot", "junction_pwm") %in% names(result)))
    expect_s3_class(result$plot, "gg")
    expect_true(is.matrix(result$junction_pwm))

    seqs <- c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF")
    p <- plot_cdr3_logo(seqs, chain = "alpha",
                         return_junction_pwm = FALSE)
    expect_s3_class(p, "gg")
})


# ===========================================================================
# plot_tcr_logo_panel
# ===========================================================================

test_that("plot_tcr_logo_panel works with various options", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")
    skip_if_not_installed("patchwork")
    skip_if_not_installed("png")
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:10, ]

    # With junction bars
    p1 <- plot_tcr_logo_panel(pa, organism = "mouse",
                               show_junction_bars = TRUE)
    expect_true(inherits(p1, "patchwork") || inherits(p1, "gg"))

    # Without junction bars
    p2 <- plot_tcr_logo_panel(pa, organism = "mouse",
                               show_junction_bars = FALSE)
    expect_true(inherits(p2, "patchwork") || inherits(p2, "gg"))

    # With title
    p3 <- plot_tcr_logo_panel(pa, organism = "mouse",
                               title = "PA panel",
                               show_junction_bars = FALSE)
    expect_true(inherits(p3, "patchwork") || inherits(p3, "gg"))

    # Missing nucseq: should still work, just no junction bars
    pa_no_nuc <- pa
    pa_no_nuc$cdr3a_nucseq <- NULL
    pa_no_nuc$cdr3b_nucseq <- NULL
    p4 <- plot_tcr_logo_panel(pa_no_nuc, organism = "mouse",
                               show_junction_bars = TRUE)
    expect_true(inherits(p4, "patchwork") || inherits(p4, "gg"))
})


# ===========================================================================
# TCRrep input support
# ===========================================================================

test_that("plot functions accept TCRrep input", {
    skip_on_cran()
    data(dash, envir = environment())
    rep <- TCRrep(dash[1:30, ], "mouse", compute_distances = TRUE)

    # Distance-based plots
    expect_true(inherits(plot_tcrdist_heatmap(tcr_rep = rep), "ggplot"))
    expect_true(inherits(plot_distance_distribution(tcr_rep = rep), "ggplot"))
    expect_true(inherits(plot_tcrdist_dendrogram(tcr_rep = rep), "ggplot"))

    # tcr_df-based plots
    expect_true(inherits(plot_gene_usage(tcr_rep = rep, gene_col = "va"),
                         "ggplot"))
    expect_true(inherits(plot_cdr3_length(tcr_rep = rep), "ggplot"))

    # Scatter with metadata from TCRrep
    coords <- matrix(rnorm(60), ncol = 2)
    expect_true(inherits(
        plot_tcr_scatter(coords, color_by = "epitope", tcr_rep = rep),
        "ggplot"))

    # Alluvial (sankey mode, no ggalluvial needed)
    expect_true(inherits(plot_gene_alluvial(tcr_rep = rep), "ggplot"))

    # Conflict errors
    mat <- rep@paired_dist
    expect_error(plot_tcrdist_heatmap(mat, tcr_rep = rep), "not both")
    expect_error(plot_distance_distribution(mat, tcr_rep = rep), "not both")
    expect_error(plot_tcrdist_dendrogram(rep@clone_df, tcr_rep = rep),
                 "not both")
    expect_error(plot_gene_usage(rep@clone_df, "va", tcr_rep = rep),
                 "not both")
    expect_error(plot_cdr3_length(rep@clone_df, tcr_rep = rep), "not both")
    expect_error(plot_gene_alluvial(rep@clone_df, tcr_rep = rep), "not both")
})


# ===========================================================================
# plot_gene_alluvial
# ===========================================================================

test_that("plot_gene_alluvial sankey mode works with all options", {
    data(dash, envir = environment())

    # Default sankey 4-axis
    p <- plot_gene_alluvial(dash)
    expect_s3_class(p, "ggplot")

    # Custom 2-axis
    p2 <- plot_gene_alluvial(dash, axes = c("va", "vb"))
    expect_s3_class(p2, "ggplot")

    # 3-axis
    p3 <- plot_gene_alluvial(dash, axes = c("ja", "va", "vb"))
    expect_s3_class(p3, "ggplot")

    # max_genes caps rare genes
    p4 <- plot_gene_alluvial(dash, max_genes = 3L)
    expect_s3_class(p4, "ggplot")

    # Facet by column name
    p5 <- plot_gene_alluvial(dash, facet_by = "epitope")
    expect_s3_class(p5, "ggplot")
    expect_false(inherits(p5$facet, "FacetNull"))

    # Facet by vector
    p6 <- plot_gene_alluvial(dash, facet_by = dash$epitope)
    expect_s3_class(p6, "ggplot")

    # Colored by axis
    p7 <- plot_gene_alluvial(dash, color_by = "va")
    expect_s3_class(p7, "ggplot")

    # Custom fill color and alpha
    p8 <- plot_gene_alluvial(dash, fill_color = "steelblue", alpha = 0.8)
    expect_s3_class(p8, "ggplot")
})

test_that("plot_gene_alluvial alluvial mode works", {
    skip_if_not_installed("ggalluvial")
    data(dash, envir = environment())

    p <- plot_gene_alluvial(dash, mode = "alluvial")
    expect_s3_class(p, "ggplot")

    p2 <- plot_gene_alluvial(dash, mode = "alluvial", color_by = "va")
    expect_s3_class(p2, "ggplot")

    p3 <- plot_gene_alluvial(dash, mode = "alluvial",
                              fill_color = "steelblue", alpha = 0.8)
    expect_s3_class(p3, "ggplot")

    p4 <- plot_gene_alluvial(dash, mode = "alluvial", facet_by = "epitope")
    expect_s3_class(p4, "ggplot")
})

test_that("plot_gene_alluvial input validation", {
    data(dash, envir = environment())

    # Invalid color_by
    expect_error(plot_gene_alluvial(dash, color_by = "epitope"),
                 "must be one of the axes")

    # Missing column error
    expect_error(plot_gene_alluvial(dash, axes = c("va", "nonexistent")),
                 "missing columns")
})
