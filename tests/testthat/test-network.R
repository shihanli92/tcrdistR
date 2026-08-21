# Tests for TCR network computation and plotting.
# Internal helpers (.find_distance_valley, .gmm_2_component, .jitter_overlapping)
# are tested implicitly through compute_tcr_network.


# ===========================================================================
# compute_tcr_network
# ===========================================================================

test_that("compute_tcr_network returns correct structure with metadata", {
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
    expect_equal(igraph::vcount(net$graph), 100L)
    # Metadata attached as vertex attributes
    vattr <- igraph::vertex_attr_names(net$graph)
    expect_true("va" %in% vattr)
    expect_true("cdr3b" %in% vattr)
    expect_true("epitope" %in% vattr)
    expect_equal(igraph::vertex_attr(net$graph, "epitope"), sub$epitope)
    # Edges have weight and distance
    if (net$n_edges > 0L) {
        expect_true("weight" %in% igraph::edge_attr_names(net$graph))
        expect_true("distance" %in% igraph::edge_attr_names(net$graph))
        w <- igraph::E(net$graph)$weight
        expect_true(all(w > 0 & w <= 1))
        d <- igraph::E(net$graph)$distance
        expect_true(all(d <= 48))
    }
})

test_that("compute_tcr_network threshold behavior and dist_matrix input", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:50, ]

    # Higher threshold connects more pairs
    net_low  <- compute_tcr_network(sub, "mouse", threshold = 20, seed = 1)
    net_high <- compute_tcr_network(sub, "mouse", threshold = 100, seed = 1)
    expect_true(net_high$n_edges >= net_low$n_edges)

    # Auto-detect threshold
    net_auto <- suppressMessages(
        compute_tcr_network(sub, "mouse", seed = 42)
    )
    expect_true(net_auto$threshold > 0 && net_auto$threshold < 500)

    # dist_matrix input
    dm <- tcrdist_matrix(sub, "mouse")
    net_dm <- compute_tcr_network(sub, dist_matrix = dm, threshold = 48,
                                   seed = 1)
    expect_equal(igraph::vcount(net_dm$graph), 50L)
    expect_equal(net_dm$threshold, 48)

    # Auto-detect from dist_matrix
    net_dm_auto <- suppressMessages(
        compute_tcr_network(sub, dist_matrix = dm, seed = 1)
    )
    expect_true(net_dm_auto$threshold > 0)
})

test_that("compute_tcr_network reproducibility and jitter", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:50, ]

    n1 <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42)
    n2 <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42)
    expect_equal(n1$layout, n2$layout)

    n3 <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42,
                               jitter = FALSE)
    n4 <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 42,
                               jitter = FALSE)
    expect_equal(n3$layout, n4$layout)
})

test_that("compute_tcr_network min_edges pruning", {
    skip_if_not_installed("igraph")
    data(dash)
    sub <- dash[1:100, ]

    net_full   <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)
    net_pruned <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1,
                                       min_edges = 1L)
    deg <- igraph::degree(net_pruned$graph)
    expect_true(all(deg >= 1L))
    expect_true(igraph::vcount(net_pruned$graph) <=
                    igraph::vcount(net_full$graph))

    # min_edges=0 keeps all vertices
    net_all <- compute_tcr_network(dash[1:30, ], "mouse", threshold = 48,
                                    seed = 1, min_edges = 0L)
    expect_equal(igraph::vcount(net_all$graph), 30L)

    # Higher min_edges threshold
    net_high <- compute_tcr_network(sub, "mouse", threshold = 100, seed = 1,
                                     min_edges = 3L)
    expect_true(all(igraph::degree(net_high$graph) >= 3L))

    # Warns if all removed
    expect_warning(
        net_warn <- compute_tcr_network(dash[1:10, ], "mouse", threshold = 1,
                                         seed = 1, min_edges = 100L),
        "removed all vertices"
    )
    expect_equal(igraph::vcount(net_warn$graph), 10L)
})

test_that("compute_tcr_network input validation", {
    skip_if_not_installed("igraph")
    data(dash)

    expect_error(compute_tcr_network(dash[1:10, ], threshold = 48),
                 "organism.*required")
    expect_error(compute_tcr_network("not_df"), "must be a data.frame")
    expect_error(compute_tcr_network(dash[1:10, ], "mouse", threshold = 48,
                                      layout = "bad"),
                 "Unknown layout")
})


# ===========================================================================
# plot_tcr_network
# ===========================================================================

test_that("plot_tcr_network handles various color_by modes", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggplot2")
    data(dash)
    sub <- dash[1:50, ]
    net <- compute_tcr_network(sub, "mouse", threshold = 48, seed = 1)

    # String attribute, vector, NULL, continuous, custom palette
    expect_true(inherits(plot_tcr_network(net, color_by = "epitope"), "ggplot"))
    expect_true(inherits(plot_tcr_network(net, color_by = sub$epitope),
                         "ggplot"))
    expect_true(inherits(plot_tcr_network(net), "ggplot"))
    expect_true(inherits(plot_tcr_network(net, color_by = seq_len(50)),
                         "ggplot"))
    expect_true(inherits(
        plot_tcr_network(net, color_by = "epitope",
                          palette = c("red", "blue", "green")), "ggplot"))
})

test_that("plot_tcr_network works with no edges", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggplot2")
    data(dash)
    net <- compute_tcr_network(dash[1:10, ], "mouse", threshold = 1, seed = 1)
    expect_true(inherits(plot_tcr_network(net, color_by = "epitope"), "ggplot"))
})

test_that("compute_tcr_network accepts TCRrep input", {
    skip_if_not_installed("igraph")
    data(dash)
    rep <- TCRrep(dash[1:30, ], "mouse", compute_distances = TRUE)

    net <- compute_tcr_network(tcr_rep = rep, threshold = 100, seed = 1)
    expect_true(inherits(net$graph, "igraph"))
    expect_equal(igraph::vcount(net$graph), 30L)

    # Conflict error
    expect_error(compute_tcr_network(rep@clone_df, tcr_rep = rep), "not both")
})
