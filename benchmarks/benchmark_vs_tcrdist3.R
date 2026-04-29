#!/usr/bin/env Rscript
#
# Benchmark tcrdistR vs Python tcrdist3
#
# Usage: Rscript benchmarks/benchmark_vs_tcrdist3.R [path_to_dash.csv]
#
# If no path is given, looks for tests/testthat/fixtures/dash.csv in the
# package source tree (assumes running from the package root).

library(tcrdistR)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) {
    dash_file <- args[1L]
} else {
    dash_file <- "tests/testthat/fixtures/dash.csv"
}

if (!file.exists(dash_file)) {
    stop("Cannot find dash.csv at: ", dash_file,
         "\nRun from the package root or provide the path as an argument.")
}

cat("Loading data from:", dash_file, "\n")
dash <- utils::read.csv(dash_file, stringsAsFactors = FALSE)
dash <- data.frame(
    va    = dash$v_a_gene,
    cdr3a = dash$cdr3_a_aa,
    vb    = dash$v_b_gene,
    cdr3b = dash$cdr3_b_aa,
    stringsAsFactors = FALSE
)
cat("Total TCRs:", nrow(dash), "\n\n")

# Warm up (gene database loading, JIT)
invisible(tcrdist_matrix(dash[1:10, ], "mouse"))

sizes <- c(50L, 100L, 200L, 500L, 1000L, nrow(dash))
sizes <- sizes[sizes <= nrow(dash)]

results <- data.frame(
    n = integer(0),
    pairs = numeric(0),
    time_sec = numeric(0),
    pairs_per_sec = numeric(0)
)

for (n in sizes) {
    sub <- dash[seq_len(n), ]
    npairs <- as.double(n) * (n - 1) / 2

    timing <- system.time({
        mat <- tcrdist_matrix(sub, "mouse")
    })

    elapsed <- timing["elapsed"]
    rate <- npairs / elapsed

    results <- rbind(results, data.frame(
        n = n,
        pairs = npairs,
        time_sec = round(elapsed, 4),
        pairs_per_sec = round(rate, 0)
    ))

    cat(sprintf("  n=%5d | %10.0f pairs | %7.4f sec | %10.0f pairs/sec\n",
                n, npairs, elapsed, rate))
}

cat("\n--- tcrdistR Benchmark Results ---\n")
print(results, row.names = FALSE)
cat("\nR version:", R.version.string, "\n")
cat("tcrdistR version:", as.character(packageVersion("tcrdistR")), "\n")
