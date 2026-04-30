# ===========================================================================
# .find_distance_valley tests
# ===========================================================================

test_that(".find_distance_valley finds valley in bimodal distribution", {
    set.seed(42)
    # Two well-separated normal peaks at 10 and 100
    dists <- c(rnorm(500, mean = 10, sd = 5), rnorm(500, mean = 100, sd = 15))
    valley <- tcrdistR:::.find_distance_valley(dists)
    # Valley should be between the peaks
    expect_true(valley > 20 && valley < 80)
})

test_that(".find_distance_valley finds valley in widely separated peaks", {
    set.seed(123)
    # Mimic real TCRdist: sharp peak near 0, broad peak at ~200
    dists <- c(rnorm(800, mean = 5, sd = 3),
               rnorm(200, mean = 200, sd = 40))
    dists <- dists[dists >= 0]  # distances are non-negative
    valley <- tcrdistR:::.find_distance_valley(dists)
    # Valley should be between the two peaks, not near the first peak
    expect_true(valley > 30 && valley < 180)
})

test_that(".find_distance_valley handles unimodal with fallback", {
    set.seed(1)
    dists <- rnorm(200, mean = 50, sd = 5)
    # GMM or KDE may handle unimodal — just check result is reasonable
    valley <- suppressMessages(
        tcrdistR:::.find_distance_valley(dists)
    )
    expect_true(is.numeric(valley))
    expect_true(valley > 0)
})

test_that(".gmm_2_component finds crossover between Gaussian components", {
    set.seed(99)
    x <- c(rnorm(500, mean = 20, sd = 8), rnorm(500, mean = 150, sd = 20))
    crossover <- tcrdistR:::.gmm_2_component(x)
    expect_true(is.numeric(crossover))
    # Should be between the two means
    expect_true(crossover > 40 && crossover < 130)
})

test_that(".gmm_2_component returns NULL for too few points", {
    expect_null(tcrdistR:::.gmm_2_component(1:10))
})

test_that(".find_distance_valley warns on too few distances", {
    expect_warning(
        valley <- tcrdistR:::.find_distance_valley(c(1, 2, 3)),
        "too few"
    )
    expect_equal(valley, 2)
})


# ===========================================================================
# .jitter_overlapping tests
# ===========================================================================

test_that(".jitter_overlapping spreads overlapping nodes", {
    coords <- matrix(c(1, 1, 1, 2,
                        3, 3, 3, 4), ncol = 2)
    # Rows 1-3 overlap at (1,3), row 4 is unique
    result <- tcrdistR:::.jitter_overlapping(coords)
    # Overlapping nodes should now have distinct coordinates
    expect_false(all(result[1, ] == result[2, ]))
    expect_false(all(result[1, ] == result[3, ]))
    # Non-overlapping node should be unchanged
    expect_equal(result[4, ], coords[4, ])
})

test_that(".jitter_overlapping is a no-op without overlaps", {
    coords <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
    result <- tcrdistR:::.jitter_overlapping(coords)
    expect_equal(result, coords)
})


# ===========================================================================
# compute_tcr_network tests
# ===========================================================================

test_that("compute_tcr_network returns correct structure", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:100, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42)

    expect_true(is.list(net))
    expect_true(igraph::is_igraph(net$graph))
    expect_true(is.matrix(net$layout))
    expect_equal(nrow(net$layout), 100L)
    expect_equal(ncol(net$layout), 2L)
    expect_equal(net$threshold, 48)
    expect_true(is.numeric(net$scale))
    expect_true(is.numeric(net$n_components))
    expect_true(is.numeric(net$n_edges))
    expect_true(inherits(net$dist_plot, "ggplot"))
})

test_that("compute_tcr_network vertex count matches input", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:50, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)
    expect_equal(igraph::vcount(net$graph), 50L)
})

test_that("compute_tcr_network attaches metadata as vertex attributes", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:30, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)

    # All columns from tcrs should be vertex attributes
    vattr <- igraph::vertex_attr_names(net$graph)
    expect_true("va" %in% vattr)
    expect_true("cdr3b" %in% vattr)
    expect_true("epitope" %in% vattr)
    expect_equal(igraph::vertex_attr(net$graph, "epitope"), sub$epitope)
})

test_that("compute_tcr_network edges have weight and distance", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:50, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)

    if (net$n_edges > 0L) {
        expect_true("weight" %in% igraph::edge_attr_names(net$graph))
        expect_true("distance" %in% igraph::edge_attr_names(net$graph))
        # Weights should be in (0, 1]
        w <- igraph::E(net$graph)$weight
        expect_true(all(w > 0 & w <= 1))
        # Distances should be <= threshold
        d <- igraph::E(net$graph)$distance
        expect_true(all(d <= 48))
    }
})

test_that("compute_tcr_network with dist_matrix input", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:30, ]

    dm <- tcrdist_matrix(sub, "mouse")
    net <- compute_tcr_network(sub, dist_matrix = dm, threshold = 48, seed = 1)

    expect_equal(igraph::vcount(net$graph), 30L)
    expect_equal(net$threshold, 48)
})

test_that("compute_tcr_network auto-detects threshold", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:100, ]

    net <- suppressMessages(
        compute_tcr_network(sub, "mouse", seed = 42)
    )
    expect_true(net$threshold > 0)
    expect_true(net$threshold < 500)
})

test_that("compute_tcr_network auto-detects threshold from dist_matrix", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:50, ]

    dm <- tcrdist_matrix(sub, "mouse")
    net <- suppressMessages(
        compute_tcr_network(sub, dist_matrix = dm, seed = 1)
    )
    expect_true(net$threshold > 0)
})

test_that("compute_tcr_network is reproducible with seed", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:50, ]

    n1 <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42)
    n2 <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42)
    expect_equal(n1$layout, n2$layout)
})

test_that("compute_tcr_network jitter=FALSE keeps raw layout", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:50, ]

    n1 <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42,
                               jitter = FALSE)
    n2 <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42,
                               jitter = FALSE)
    expect_equal(n1$layout, n2$layout)
})

test_that("compute_tcr_network errors without organism or dist_matrix", {
    skip_if_not_installed("igraph")
    data(dash)
    expect_error(
        compute_tcr_network(dash[1:10, ], threshold = 48),
        "organism.*required"
    )
})

test_that("compute_tcr_network errors on non-data.frame", {
    skip_if_not_installed("igraph")
    expect_error(compute_tcr_network("not_df"), "must be a data.frame")
})

test_that("compute_tcr_network errors on bad layout", {
    skip_if_not_installed("igraph")
    data(dash)
    expect_error(
        compute_tcr_network(dash[1:10, ], "mouse", threshold = 48,
                             layout = "bad"),
        "Unknown layout"
    )
})

test_that("compute_tcr_network high threshold connects more", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:50, ]

    net_low  <- compute_tcr_network(sub, "mouse", threshold = 20, seed = 1)
    net_high <- compute_tcr_network(sub, "mouse", threshold = 100, seed = 1)

    expect_true(net_high$n_edges >= net_low$n_edges)
})

test_that("compute_tcr_network min_edges=1 removes singletons", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:50, ]

    net_full   <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)
    net_pruned <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1,
                                       min_edges = 1L)

    # Pruned graph should have no isolated nodes
    deg <- igraph::degree(net_pruned$graph)
    expect_true(all(deg >= 1L))
    # And fewer or equal vertices
    expect_true(igraph::vcount(net_pruned$graph) <=
                    igraph::vcount(net_full$graph))
})

test_that("compute_tcr_network min_edges=0 keeps all vertices", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:30, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1,
                                min_edges = 0L)
    expect_equal(igraph::vcount(net$graph), 30L)
})

test_that("compute_tcr_network min_edges prunes to higher degree", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:100, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 100, seed = 1,
                                min_edges = 3L)
    deg <- igraph::degree(net$graph)
    expect_true(all(deg >= 3L))
})

test_that("compute_tcr_network min_edges warns if all removed", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:10, ]

    # Very low threshold + high min_edges = likely removes everything
    expect_warning(
        net <- compute_tcr_network(sub, "mouse", threshold = 1, seed = 1,
                                    min_edges = 100L),
        "removed all vertices"
    )
    # Should return unpruned graph
    expect_equal(igraph::vcount(net$graph), 10L)
})


# ===========================================================================
# plot_tcr_network tests
# ===========================================================================

test_that("plot_tcr_network returns ggplot with color_by string", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggplot2")
    data(dash)
    sub <- dash[1:50, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)
    p <- plot_tcr_network(net, color_by = "epitope")
    expect_true(inherits(p, "ggplot"))
})

test_that("plot_tcr_network returns ggplot with color_by vector", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggplot2")
    data(dash)
    sub <- dash[1:50, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)
    p <- plot_tcr_network(net, color_by = sub$epitope)
    expect_true(inherits(p, "ggplot"))
})

test_that("plot_tcr_network returns ggplot with NULL color_by", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggplot2")
    data(dash)
    sub <- dash[1:50, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)
    p <- plot_tcr_network(net)
    expect_true(inherits(p, "ggplot"))
})

test_that("plot_tcr_network works with continuous color", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggplot2")
    data(dash)
    sub <- dash[1:50, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)
    p <- plot_tcr_network(net, color_by = seq_len(50))
    expect_true(inherits(p, "ggplot"))
})

test_that("plot_tcr_network custom palette", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggplot2")
    data(dash)
    sub <- dash[1:50, ]

    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)
    p <- plot_tcr_network(net, color_by = "epitope",
                           palette = c("red", "blue", "green"))
    expect_true(inherits(p, "ggplot"))
})

test_that("plot_tcr_network with no edges still works", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggplot2")
    data(dash)
    sub <- dash[1:10, ]

    # Very low threshold = likely no edges
    net <- compute_tcr_network(sub, "mouse", threshold = 1, seed = 1)
    p <- plot_tcr_network(net, color_by = "epitope")
    expect_true(inherits(p, "ggplot"))
})
