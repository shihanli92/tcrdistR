# Tests for Phase 6: Visualization functions
# Tests verify structure and correctness, not visual appearance.


# ===========================================================================
# BLOSUM62 matrix
# ===========================================================================

test_that("BLOSUM62 matrix has correct dimensions and is symmetric", {
    expect_true(exists("BLOSUM62", envir = asNamespace("tcrdistR")))
    mat <- tcrdistR:::BLOSUM62
    expect_equal(nrow(mat), 20L)
    expect_equal(ncol(mat), 20L)
    expect_identical(rownames(mat), AMINO_ACIDS)
    expect_identical(colnames(mat), AMINO_ACIDS)
    expect_true(isSymmetric(mat))
})

test_that("BLOSUM62 diagonal values match known reference", {
    mat <- tcrdistR:::BLOSUM62
    expect_equal(mat["A", "A"], 4L)
    expect_equal(mat["W", "W"], 11L)
    expect_equal(mat["C", "C"], 9L)
    expect_equal(mat["Y", "Y"], 7L)
})

test_that("BLOSUM62 off-diagonal values match known reference", {
    mat <- tcrdistR:::BLOSUM62
    expect_equal(mat["A", "C"], 0L)
    expect_equal(mat["C", "A"], 0L)
    expect_equal(mat["D", "E"], 2L)
    expect_equal(mat["W", "Y"], 2L)
    expect_equal(mat["F", "Y"], 3L)
})


# ===========================================================================
# Internal: .tcrdistR_palette
# ===========================================================================

test_that(".tcrdistR_palette returns correct number of colors", {
    pal5  <- tcrdistR:::.tcrdistR_palette(5)
    pal10 <- tcrdistR:::.tcrdistR_palette(10)
    pal15 <- tcrdistR:::.tcrdistR_palette(15)
    expect_length(pal5, 5L)
    expect_length(pal10, 10L)
    expect_length(pal15, 15L)
})

test_that(".tcrdistR_palette returns valid hex colors", {
    pal <- tcrdistR:::.tcrdistR_palette(20)
    expect_true(all(grepl("^#[0-9a-fA-F]{6}$", pal)))
})

test_that(".tcrdistR_palette uses tab10 for n<=10, tab20 for n>10", {
    pal10 <- tcrdistR:::.tcrdistR_palette(10)
    pal11 <- tcrdistR:::.tcrdistR_palette(11)
    # tab10 first color is "#1f77b4"
    expect_equal(pal10[1L], "#1f77b4")
    # tab20 second color is "#aec7e8" (lighter variant)
    expect_equal(pal11[2L], "#aec7e8")
})


# ===========================================================================
# Internal: .align_cdr3_regions
# ===========================================================================

test_that(".align_cdr3_regions returns unchanged for equal length", {
    result <- tcrdistR:::.align_cdr3_regions("ABCDE", "FGHIJ")
    expect_equal(result$a, "ABCDE")
    expect_equal(result$b, "FGHIJ")
    expect_equal(nchar(result$a), nchar(result$b))
})

test_that(".align_cdr3_regions inserts gaps in shorter sequence", {
    result <- tcrdistR:::.align_cdr3_regions("CASSI", "CASSILY")
    expect_equal(nchar(result$a), nchar(result$b))
    # The shorter one should have gaps
    expect_true(grepl("-", result$a))
    expect_false(grepl("-", result$b))
    # Gap count equals length difference
    n_gaps <- nchar(gsub("[^-]", "", result$a))
    expect_equal(n_gaps, nchar("CASSILY") - nchar("CASSI"))
})

test_that(".align_cdr3_regions handles single-char shorter", {
    result <- tcrdistR:::.align_cdr3_regions("A", "ABC")
    expect_equal(nchar(result$a), 3L)
    expect_equal(nchar(result$b), 3L)
})

test_that(".align_cdr3_regions custom gap character", {
    result <- tcrdistR:::.align_cdr3_regions("CASSI", "CASSILY",
                                              gap_character = ".")
    expect_true(grepl("\\.", result$a))
})


# ===========================================================================
# Internal: .build_cdr3_pwm
# ===========================================================================

test_that(".build_cdr3_pwm single sequence gives identity PWM", {
    result <- tcrdistR:::.build_cdr3_pwm("CASSIRSSYEQY", trim = FALSE)
    expect_true(is.matrix(result$pwm))
    expect_equal(result$center_idx, 1L)
    # Each column should have exactly one 1.0
    for (j in seq_len(ncol(result$pwm))) {
        expect_equal(sum(result$pwm[, j]), 1.0, tolerance = 1e-10)
        expect_equal(max(result$pwm[, j]), 1.0)
    }
})

test_that(".build_cdr3_pwm columns sum to 1", {
    seqs <- c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSSYEQY",
              "CASSIRSSYEQF")
    result <- tcrdistR:::.build_cdr3_pwm(seqs, trim = FALSE)
    col_sums <- colSums(result$pwm)
    expect_true(all(abs(col_sums - 1.0) < 1e-10))
})

test_that(".build_cdr3_pwm center_idx is valid", {
    seqs <- c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY")
    result <- tcrdistR:::.build_cdr3_pwm(seqs, trim = FALSE)
    expect_true(result$center_idx >= 1L)
    expect_true(result$center_idx <= length(seqs))
})

test_that(".build_cdr3_pwm trim removes 3+2 residues", {
    seq_full <- "CASSIRSSYEQYF"  # 13 chars -> trim to 8 (13-3-2)
    result_trim <- tcrdistR:::.build_cdr3_pwm(seq_full, trim = TRUE)
    result_notrim <- tcrdistR:::.build_cdr3_pwm(seq_full, trim = FALSE)
    expect_equal(ncol(result_trim$pwm), nchar(seq_full) - 5L)
    expect_equal(ncol(result_notrim$pwm), nchar(seq_full))
})

test_that(".build_cdr3_pwm aligned_seqs have equal length", {
    seqs <- c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY")
    result <- tcrdistR:::.build_cdr3_pwm(seqs, trim = FALSE)
    lengths <- nchar(result$aligned_seqs)
    expect_true(length(unique(lengths)) == 1L)
})


# ===========================================================================
# Internal: .hclust_to_segments
# ===========================================================================

test_that(".hclust_to_segments returns correct structure", {
    set.seed(42)
    d <- stats::dist(matrix(rnorm(20), ncol = 2))
    hc <- stats::hclust(d)
    seg <- tcrdistR:::.hclust_to_segments(hc)
    expect_s3_class(seg, "data.frame")
    expect_true(all(c("x", "y", "xend", "yend") %in% colnames(seg)))
})

test_that(".hclust_to_segments all y values >= 0", {
    set.seed(42)
    d <- stats::dist(matrix(rnorm(20), ncol = 2))
    hc <- stats::hclust(d)
    seg <- tcrdistR:::.hclust_to_segments(hc)
    expect_true(all(seg$y >= 0))
    expect_true(all(seg$yend >= 0))
})

test_that(".hclust_to_segments correct segment count", {
    # For n leaves: n-1 horizontal + 2*(n-1) vertical = 3*(n-1)
    set.seed(42)
    n <- 5L
    d <- stats::dist(matrix(rnorm(n * 2), ncol = 2))
    hc <- stats::hclust(d)
    seg <- tcrdistR:::.hclust_to_segments(hc)
    # Each merge generates 2 vertical + 1 horizontal = 3 segments
    expect_equal(nrow(seg), 3L * (n - 1L))
})


# ===========================================================================
# Exported: plot functions (require ggplot2)
# ===========================================================================

test_that("plot_cdr3_length returns ggplot for alpha chain", {
    skip_if_not_installed("ggplot2")

    tcr_df <- data.frame(
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF",
                   "CAVRDSSYKLIF", "CAVKDSYKLIF"),
        stringsAsFactors = FALSE
    )
    p <- plot_cdr3_length(tcr_df, chain = "alpha")
    expect_s3_class(p, "gg")
})

test_that("plot_cdr3_length returns ggplot for both chains", {
    skip_if_not_installed("ggplot2")

    tcr_df <- data.frame(
        cdr3a = c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF"),
        cdr3b = c("CASSIRSSYEQY", "CASSIKSSYEQY", "CASSIRSYEQY"),
        stringsAsFactors = FALSE
    )
    p <- plot_cdr3_length(tcr_df, chain = "both")
    expect_s3_class(p, "gg")
})

test_that("plot_cdr3_length handles empty CDR3 gracefully", {
    skip_if_not_installed("ggplot2")

    tcr_df <- data.frame(
        cdr3a = character(0),
        cdr3b = character(0),
        stringsAsFactors = FALSE
    )
    p <- plot_cdr3_length(tcr_df, chain = "both")
    expect_s3_class(p, "gg")
})


test_that("plot_gene_usage returns ggplot", {
    skip_if_not_installed("ggplot2")

    tcr_df <- data.frame(
        va = c("TRAV1-2*01", "TRAV1-2*02", "TRAV12-1*01",
               "TRAV1-2*01", "TRAV12-1*01"),
        stringsAsFactors = FALSE
    )
    p <- plot_gene_usage(tcr_df, "va")
    expect_s3_class(p, "gg")
})

test_that("plot_gene_usage handles empty data", {
    skip_if_not_installed("ggplot2")

    tcr_df <- data.frame(va = character(0), stringsAsFactors = FALSE)
    p <- plot_gene_usage(tcr_df, "va")
    expect_s3_class(p, "gg")
})

test_that("plot_gene_usage max_genes limits bars", {
    skip_if_not_installed("ggplot2")

    genes <- paste0("TRAV", seq_len(30))
    tcr_df <- data.frame(va = genes, stringsAsFactors = FALSE)
    p <- plot_gene_usage(tcr_df, "va", strip_allele = FALSE, max_genes = 5L)
    # Build the plot to extract data
    pb <- ggplot2::ggplot_build(p)
    expect_true(nrow(pb$data[[1L]]) <= 5L)
})

test_that("plot_tcrdist_heatmap returns ggplot", {
    skip_if_not_installed("ggplot2")

    mat <- matrix(c(0, 10, 20,
                     10, 0, 15,
                     20, 15, 0), nrow = 3)
    p <- plot_tcrdist_heatmap(mat)
    expect_s3_class(p, "gg")
})

test_that("plot_tcrdist_heatmap cluster=FALSE works", {
    skip_if_not_installed("ggplot2")

    mat <- matrix(c(0, 10, 20,
                     10, 0, 15,
                     20, 15, 0), nrow = 3)
    p <- plot_tcrdist_heatmap(mat, cluster = FALSE,
                               labels = c("A", "B", "C"))
    expect_s3_class(p, "gg")
})

test_that("plot_distance_distribution returns ggplot", {
    skip_if_not_installed("ggplot2")

    mat <- matrix(c(0, 10, 20, 30,
                     10, 0, 15, 25,
                     20, 15, 0, 12,
                     30, 25, 12, 0), nrow = 4)
    p <- plot_distance_distribution(mat)
    expect_s3_class(p, "gg")
})

test_that("plot_distance_distribution with custom binwidth", {
    skip_if_not_installed("ggplot2")

    mat <- matrix(c(0, 10, 20, 30,
                     10, 0, 15, 25,
                     20, 15, 0, 12,
                     30, 25, 12, 0), nrow = 4)
    p <- plot_distance_distribution(mat, binwidth = 5)
    expect_s3_class(p, "gg")
})

test_that("plot_distance_distribution with threshold line", {
    skip_if_not_installed("ggplot2")

    mat <- matrix(c(0, 10, 20, 30,
                     10, 0, 15, 25,
                     20, 15, 0, 12,
                     30, 25, 12, 0), nrow = 4)
    p <- plot_distance_distribution(mat, threshold = 18)
    expect_s3_class(p, "gg")
    # Should have the annotation layer for the threshold label
    expect_true(length(p$layers) > 3L)
})

test_that("plot_distance_distribution accepts a numeric vector", {
    skip_if_not_installed("ggplot2")

    dists <- c(5, 10, 15, 20, 25, 30)
    p <- plot_distance_distribution(dists, threshold = 12)
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter returns ggplot with NULL color_by", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    p <- plot_tcr_scatter(coords)
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter continuous coloring", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    color_by <- runif(10)
    p <- plot_tcr_scatter(coords, color_by = color_by)
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter categorical coloring", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    color_by <- rep(c("A", "B"), 5)
    p <- plot_tcr_scatter(coords, color_by = color_by)
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter custom palette", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    color_by <- rep(c("A", "B"), 5)
    p <- plot_tcr_scatter(coords, color_by = color_by,
                           palette = c("red", "blue"))
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter facet_by with vector", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    p <- plot_tcr_scatter(coords, facet_by = rep(c("A", "B"), 5))
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter facet_by with 2-column data.frame", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    facets <- data.frame(row = rep(c("R1", "R2"), 5),
                          col = rep(c("C1", "C2"), each = 5))
    p <- plot_tcr_scatter(coords, color_by = rep(c("X", "Y"), 5),
                           facet_by = facets)
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter highlight with logical vector", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    hl <- c(TRUE, TRUE, TRUE, rep(FALSE, 7))
    p <- plot_tcr_scatter(coords, highlight = hl)
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter highlight with integer indices", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    p <- plot_tcr_scatter(coords, highlight = c(1L, 5L, 10L))
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter highlight with color_by", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    hl <- c(TRUE, TRUE, TRUE, rep(FALSE, 7))
    p <- plot_tcr_scatter(coords, color_by = rep(c("A", "B"), 5),
                           highlight = hl)
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter metadata column-name lookup", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    meta <- data.frame(group = rep(c("A", "B"), 5),
                        stringsAsFactors = FALSE)
    p <- plot_tcr_scatter(coords, color_by = "group", facet_by = "group",
                           metadata = meta)
    expect_s3_class(p, "gg")
})

test_that("plot_tcr_scatter metadata named-list highlight", {
    skip_if_not_installed("ggplot2")

    coords <- matrix(rnorm(20), ncol = 2)
    meta <- data.frame(group = rep(c("A", "B"), 5),
                        site = rep(c("X", "Y", "X", "Y", "X"), 2),
                        stringsAsFactors = FALSE)
    p <- plot_tcr_scatter(coords, highlight = list(group = "A", site = "X"),
                           metadata = meta)
    expect_s3_class(p, "gg")
})


# ===========================================================================
# CDR3 logo tests (require ggseqlogo)
# ===========================================================================

test_that("plot_cdr3_logo returns ggplot with bits method", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")

    seqs <- c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF",
              "CAVRDSSYKLIF", "CAVKDSYKLIF")
    p <- plot_cdr3_logo(seqs, chain = "alpha", method = "bits")
    expect_s3_class(p, "gg")
})

test_that("plot_cdr3_logo returns ggplot with prob method", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")

    seqs <- c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF",
              "CAVRDSSYKLIF", "CAVKDSYKLIF")
    p <- plot_cdr3_logo(seqs, chain = "beta", method = "prob")
    expect_s3_class(p, "gg")
})

test_that("plot_cdr3_logo works with single sequence", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")

    p <- plot_cdr3_logo("CAVRDSSYKLIF", chain = "alpha")
    expect_s3_class(p, "gg")
})


# ===========================================================================
# Dendrogram (requires gene DB + compiled code)
# ===========================================================================

test_that("plot_tcrdist_dendrogram works with small synthetic TCRs", {
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

    p <- plot_tcrdist_dendrogram(tcr_df, "human")
    expect_s3_class(p, "gg")
})

test_that("plot_tcrdist_dendrogram with categorical color_by", {
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

    p <- plot_tcrdist_dendrogram(tcr_df, "human",
                                  color_by = c("A", "A", "B", "B", "A"))
    expect_s3_class(p, "gg")
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
# .render_gene_text_raster
# ===========================================================================

test_that(".render_gene_text_raster returns RGBA array", {
    skip_if_not_installed("png")
    skip_on_cran()
    img <- tcrdistR:::.render_gene_text_raster("1-2", col = "red")
    expect_true(is.array(img))
    expect_equal(length(dim(img)), 3L)
    expect_equal(dim(img)[3L], 4L)
})

test_that(".render_gene_text_raster has non-zero dimensions", {
    skip_if_not_installed("png")
    skip_on_cran()
    img <- tcrdistR:::.render_gene_text_raster("TRAV1")
    expect_true(dim(img)[1L] > 0L)
    expect_true(dim(img)[2L] > 0L)
})


# ===========================================================================
# plot_vj_gene_logo
# ===========================================================================

test_that("plot_vj_gene_logo returns ggplot for V genes", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("png")
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ]
    p <- plot_vj_gene_logo(pa$vb, organism = "mouse",
                            gene_type = "V", chain = "beta")
    expect_s3_class(p, "gg")
})

test_that("plot_vj_gene_logo returns ggplot for J genes", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("png")
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ]
    p <- plot_vj_gene_logo(pa$ja, organism = "mouse",
                            gene_type = "J", chain = "alpha")
    expect_s3_class(p, "gg")
})

test_that("plot_vj_gene_logo respects max_genes", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("png")
    skip_on_cran()
    data(dash, envir = environment())
    p <- plot_vj_gene_logo(dash$vb, organism = "mouse",
                            gene_type = "V", chain = "beta",
                            max_genes = 3L)
    expect_s3_class(p, "gg")
})

test_that("plot_vj_gene_logo works with human organism", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("png")
    skip_on_cran()
    genes <- rep(c("TRAV1-2*01", "TRAV12-1*01", "TRAV10*01"), c(5, 3, 2))
    p <- plot_vj_gene_logo(genes, organism = "human",
                            gene_type = "V", chain = "alpha")
    expect_s3_class(p, "gg")
})


# ===========================================================================
# compute_nucseq_src
# ===========================================================================

test_that("compute_nucseq_src returns list of correct length", {
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:5, ]
    src <- compute_nucseq_src(pa, organism = "mouse", chain = "alpha")
    expect_type(src, "list")
    expect_length(src, 5L)
})

test_that("compute_nucseq_src alpha chain uses V/N/J labels", {
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:5, ]
    src <- compute_nucseq_src(pa, organism = "mouse", chain = "alpha")
    for (s in src) {
        if (!is.null(s)) {
            expect_true(all(s %in% c("V", "N", "J")))
        }
    }
})

test_that("compute_nucseq_src beta chain uses N1/N2/D labels", {
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:5, ]
    src <- compute_nucseq_src(pa, organism = "mouse", chain = "beta")
    for (s in src) {
        if (!is.null(s)) {
            expect_true(all(s %in% c("V", "N1", "D", "N2", "J")))
        }
    }
})

test_that("compute_nucseq_src returns NULL for missing nucseq", {
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:3, ]
    pa$cdr3b_nucseq <- NA_character_
    src <- compute_nucseq_src(pa, organism = "mouse", chain = "beta")
    expect_true(all(vapply(src, is.null, logical(1L))))
})


# ===========================================================================
# plot_cdr3_logo with return_junction_pwm
# ===========================================================================

test_that("plot_cdr3_logo return_junction_pwm returns list", {
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
    expect_true("plot" %in% names(result))
    expect_true("junction_pwm" %in% names(result))
    expect_s3_class(result$plot, "gg")
    expect_true(is.matrix(result$junction_pwm))
})

test_that("plot_cdr3_logo return_junction_pwm FALSE is backward compatible", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")
    skip_on_cran()
    seqs <- c("CAVRDSSYKLIF", "CAVKDSSYKLIF", "CAVRDSYKLIF")
    p <- plot_cdr3_logo(seqs, chain = "alpha",
                         return_junction_pwm = FALSE)
    expect_s3_class(p, "gg")
})


# ===========================================================================
# plot_tcr_logo_panel
# ===========================================================================

test_that("plot_tcr_logo_panel returns patchwork with junction bars", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")
    skip_if_not_installed("patchwork")
    skip_if_not_installed("png")
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:10, ]
    p <- plot_tcr_logo_panel(pa, organism = "mouse",
                              show_junction_bars = TRUE)
    expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})

test_that("plot_tcr_logo_panel works without junction bars", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")
    skip_if_not_installed("patchwork")
    skip_if_not_installed("png")
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:10, ]
    p <- plot_tcr_logo_panel(pa, organism = "mouse",
                              show_junction_bars = FALSE)
    expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})

test_that("plot_tcr_logo_panel adds title", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")
    skip_if_not_installed("patchwork")
    skip_if_not_installed("png")
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:10, ]
    p <- plot_tcr_logo_panel(pa, organism = "mouse",
                              title = "PA panel",
                              show_junction_bars = FALSE)
    expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})

test_that("plot_tcr_logo_panel gracefully handles missing nucseq", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggseqlogo")
    skip_if_not_installed("patchwork")
    skip_if_not_installed("png")
    skip_on_cran()
    data(dash, envir = environment())
    pa <- dash[dash$epitope == "PA", ][1:10, ]
    # Remove nucseq columns — should still work, just no junction bars
    pa$cdr3a_nucseq <- NULL
    pa$cdr3b_nucseq <- NULL
    p <- plot_tcr_logo_panel(pa, organism = "mouse",
                              show_junction_bars = TRUE)
    expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})
