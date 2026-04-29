test_that("load_gene_database returns genes for human", {
    genes <- load_gene_database("human")
    expect_true(length(genes) > 0)
    expect_true(is.list(genes))
})

test_that("load_gene_database returns genes for mouse", {
    genes <- load_gene_database("mouse")
    expect_true(length(genes) > 0)
})

test_that("gene entries have required fields", {
    genes <- load_gene_database("human")
    g <- genes[[1]]
    expect_true(all(c("id", "organism", "chain", "region") %in% names(g)))
})

test_that("V genes have CDRs", {
    genes <- load_gene_database("human")
    v_genes <- Filter(function(g) g$region == "V", genes)
    expect_true(length(v_genes) > 0)
    g <- v_genes[[1]]
    expect_true("cdrs" %in% names(g))
    expect_true(length(g$cdrs) > 0)
})

test_that("load_gene_database caches results", {
    # First call loads, second should be cached
    genes1 <- load_gene_database("human")
    genes2 <- load_gene_database("human")
    expect_identical(genes1, genes2)
})

test_that("trim_allele_to_gene works correctly", {
    expect_equal(trim_allele_to_gene("TRAV1-1*01"), "TRAV1-1")
    expect_equal(trim_allele_to_gene("TRBV19*01"), "TRBV19")
    expect_equal(trim_allele_to_gene("TRAV1-1"), "TRAV1-1")
})

test_that("gene database contains expected organisms", {
    human_genes <- load_gene_database("human")
    mouse_genes <- load_gene_database("mouse")
    expect_true(length(human_genes) > 0)
    expect_true(length(mouse_genes) > 0)
})
