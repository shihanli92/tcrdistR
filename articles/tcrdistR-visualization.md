# Visualization Guide

## Setup

tcrdistR provides plotting functions built on ggplot2. All return ggplot
objects that can be further customized.

``` r
library(tcrdistR)
data(dash)

# Use 100 PA-specific TCRs for fast plots
pa <- dash[dash$epitope == "PA", ]
pa_sub <- pa[1:100, ]
```

## Distance Heatmap

Visualize the pairwise distance matrix as a heatmap. Note that
[`plot_tcrdist_heatmap()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_heatmap.md)
takes a precomputed distance matrix, not a TCR data.frame:

``` r
library(ggplot2)

dist_mat <- tcrdist_matrix(pa_sub, organism = "mouse")
plot_tcrdist_heatmap(dist_mat)
```

![](tcrdistR-visualization_files/figure-html/heatmap-1.png)

With hierarchical clustering reordering (the default) and custom labels:

``` r
# Use a small subset with unique labels for readability
small_mat <- dist_mat[1:20, 1:20]
plot_tcrdist_heatmap(
  small_mat,
  labels = paste0(pa_sub$cdr3b[1:20], " (", seq_len(20), ")"),
  cluster = TRUE,
  title = "PA-specific TCR Distances"
)
```

![](tcrdistR-visualization_files/figure-html/heatmap-labels-1.png)

## Dendrogram

Plot a hierarchical clustering dendrogram. Unlike the heatmap,
[`plot_tcrdist_dendrogram()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_dendrogram.md)
takes a TCR data.frame and computes distances internally:

``` r
plot_tcrdist_dendrogram(pa_sub, organism = "mouse")
```

![](tcrdistR-visualization_files/figure-html/dendrogram-1.png)

Color by a grouping variable:

``` r
# Mix epitopes for a more interesting dendrogram
mixed <- dash[dash$epitope %in% c("PA", "NP"), ]
mixed_sub <- mixed[sample(nrow(mixed), 80), ]

plot_tcrdist_dendrogram(
  mixed_sub, organism = "mouse",
  color_by = mixed_sub$epitope,
  title = "PA vs NP TCRs"
)
```

![](tcrdistR-visualization_files/figure-html/dendrogram-colored-1.png)

## Distance Distribution

Histogram and density of pairwise distances. Like the heatmap, this
takes a precomputed distance matrix:

``` r
plot_distance_distribution(dist_mat)
```

![](tcrdistR-visualization_files/figure-html/dist-distrib-1.png)

With custom bin width:

``` r
plot_distance_distribution(dist_mat, binwidth = 10,
                           title = "PA TCR distance distribution")
```

![](tcrdistR-visualization_files/figure-html/dist-distrib-bins-1.png)

## PCA Scatter Plot

2D scatter plot from kernel PCA or any dimensionality reduction.
[`plot_tcr_scatter()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_scatter.md)
takes a N x 2 coordinate matrix:

``` r
pca <- compute_tcrdist_kernel_pca(
  mixed_sub, organism = "mouse", n_components = 5L
)

plot_tcr_scatter(
  pca$embeddings[, 1:2],
  color_by = mixed_sub$epitope,
  title = "Kernel PCA: PA vs NP",
  point_size = 2
)
```

![](tcrdistR-visualization_files/figure-html/scatter-1.png)

Color by a continuous variable (e.g., third principal component):

``` r
plot_tcr_scatter(
  pca$embeddings[, 1:2],
  color_by = pca$embeddings[, 3],
  title = "Colored by KPC3",
  point_size = 2,
  legend_title = "KPC3"
)
```

![](tcrdistR-visualization_files/figure-html/scatter-continuous-1.png)

## Gene Usage

Horizontal bar chart of V-gene frequencies.
[`plot_gene_usage()`](https://shihanli92.github.io/tcrdistR/reference/plot_gene_usage.md)
takes a data.frame and a column name:

``` r
plot_gene_usage(pa, "va", title = "V-alpha usage (PA)")
```

![](tcrdistR-visualization_files/figure-html/gene-usage-1.png)

``` r
plot_gene_usage(pa, "vb", title = "V-beta usage (PA)")
```

![](tcrdistR-visualization_files/figure-html/gene-usage-2.png)

## V/J Gene Logos

[`plot_vj_gene_logo()`](https://shihanli92.github.io/tcrdistR/reference/plot_vj_gene_logo.md)
renders gene names as text glyphs whose height is proportional to
frequency — a “sequence logo” style for gene usage:

``` r
plot_vj_gene_logo(pa$vb, organism = "mouse",
                  gene_type = "V", chain = "beta")
```

![](tcrdistR-visualization_files/figure-html/gene-logo-1.png)

``` r
plot_vj_gene_logo(pa$ja, organism = "mouse",
                  gene_type = "J", chain = "alpha")
```

![](tcrdistR-visualization_files/figure-html/gene-logo-2.png)

## CDR3 Sequence Logos

Display amino acid frequency at each position as a sequence logo. The
logo construction aligns sequences using BLOSUM62-optimal gap placement:

``` r
plot_cdr3_logo(pa_sub$cdr3b, chain = "beta", method = "bits",
               title = "CDR3-beta logo (PA)")
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the ggseqlogo package.
#>   Please report the issue at <https://github.com/omarwagih/ggseqlogo/issues>.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

![](tcrdistR-visualization_files/figure-html/logo-1.png)

Frequency (“prob”) mode:

``` r
plot_cdr3_logo(pa_sub$cdr3a, chain = "alpha", method = "prob",
               title = "CDR3-alpha logo (PA)")
```

![](tcrdistR-visualization_files/figure-html/logo-prob-1.png)

## Junction Bars

[`compute_nucseq_src()`](https://shihanli92.github.io/tcrdistR/reference/compute_nucseq_src.md)
analyzes CDR3 nucleotide sequences to identify V/N/D/J segment origins.
Pass the result to
[`plot_cdr3_logo()`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md)
with `show_junction_bars = TRUE` to display the rearrangement structure:

``` r
src <- compute_nucseq_src(pa_sub, organism = "mouse", chain = "beta")
plot_cdr3_logo(pa_sub$cdr3b, chain = "beta",
               nucseq_src = src, show_junction_bars = TRUE,
               title = "CDR3-beta with junction bars (PA)")
```

![](tcrdistR-visualization_files/figure-html/junction-1.png)

## TCR Logo Panel

[`plot_tcr_logo_panel()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_logo_panel.md)
combines V-gene logos, CDR3 sequence logos with junction bars, and
J-gene logos for both chains into a single composite panel:

``` r
plot_tcr_logo_panel(pa_sub, organism = "mouse",
                    title = "PA-specific TCR rearrangement")
```

![](tcrdistR-visualization_files/figure-html/logo-panel-1.png)

## TCR Network Plot

Visualize TCR similarity as a network graph. Nodes are TCRs, edges
connect TCRs within a distance threshold. Requires the igraph package.

``` r
# Build network with auto-detected threshold
net <- compute_tcr_network(pa_sub, organism = "mouse")
#> Auto-detected distance threshold: 216.7
```

![](tcrdistR-visualization_files/figure-html/network-1.png)

``` r

# Color by a vertex attribute (stored from the input data.frame)
plot_tcr_network(net, color_by = "subject",
                 vertex_size = 3, title = "PA TCR Network")
```

![](tcrdistR-visualization_files/figure-html/network-2.png)

Control the threshold and prune isolated nodes:

``` r
net2 <- compute_tcr_network(pa_sub, organism = "mouse",
                            threshold = 48, min_edges = 2)
```

![](tcrdistR-visualization_files/figure-html/network-custom-1.png)

``` r
plot_tcr_network(net2, vertex_size = 3,
                 title = "Threshold=48, min_edges=2")
```

![](tcrdistR-visualization_files/figure-html/network-custom-2.png)

## Combining Plots

Use the patchwork package to arrange multiple plots:

``` r
library(patchwork)

p1 <- plot_tcrdist_heatmap(dist_mat[1:30, 1:30], title = "Heatmap")
p2 <- plot_distance_distribution(dist_mat, title = "Distribution")
p3 <- plot_gene_usage(pa, "va", title = "V-alpha")
p4 <- plot_gene_usage(pa, "vb", title = "V-beta")

(p1 | p2) / (p3 | p4) +
  plot_annotation(title = "PA Repertoire Overview")
```

![](tcrdistR-visualization_files/figure-html/patchwork-1.png)
