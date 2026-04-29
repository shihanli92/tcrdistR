# Visualization Guide

## Setup

tcrdistR provides seven plotting functions built on ggplot2. All return
ggplot objects that can be further customized.

``` r
library(tcrdistR)
library(ggplot2)

tcrs <- data.frame(
  va    = c("TRAV7-3*01",  "TRAV6D-6*01", "TRAV6D-6*01",
            "TRAV6-4*01",  "TRAV6-4*01",  "TRAV7-3*01",
            "TRAV6D-6*01", "TRAV6-4*01",  "TRAV7-3*01",
            "TRAV6D-6*01"),
  cdr3a = c("CAVSLDSNYQLIW",  "CALGDRATGGNNKLTF", "CALGSNTGYQNFYF",
            "CALAPSNTNKVVF",  "CALVPSNTNKVVF",    "CAVSLDSNYQLIW",
            "CALGDRATGGNNKLTF","CALAPSNTNKVVF",    "CAVSLDSNYQLIW",
            "CALGSNTGYQNFYF"),
  vb    = c("TRBV13-1*01", "TRBV29*01",   "TRBV29*01",
            "TRBV2*01",    "TRBV29*01",   "TRBV13-1*01",
            "TRBV29*01",   "TRBV2*01",    "TRBV13-1*01",
            "TRBV29*01"),
  cdr3b = c("CASSDFDWGGDAETLYF", "CASSPDRGEVFF",    "CASTGGGAPLF",
            "CASSQDPGDYEQYF",   "CASSLGGENTLYF",   "CASSDFDWGGDAETLYF",
            "CASSPDRGEVFF",      "CASSQDPGDYEQYF",  "CASSDFDWGGDAETLYF",
            "CASTGGGAPLF"),
  epitope = c("PA", "PA", "PA", "PA", "PA",
              "NP", "NP", "NP", "NP", "NP"),
  stringsAsFactors = FALSE
)
```

## Distance Heatmap

Visualize the pairwise distance matrix as a heatmap:

``` r
# Basic heatmap
plot_tcrdist_heatmap(tcrs, organism = "mouse")

# With hierarchical clustering reordering
plot_tcrdist_heatmap(tcrs, organism = "mouse", cluster = TRUE)

# Custom color scale
plot_tcrdist_heatmap(tcrs, organism = "mouse") +
  scale_fill_viridis_c(option = "magma")
```

## Dendrogram

Plot a hierarchical clustering dendrogram:

``` r
# Basic dendrogram
plot_tcrdist_dendrogram(tcrs, organism = "mouse")

# Colored by a grouping variable
plot_tcrdist_dendrogram(tcrs, organism = "mouse",
                        color_by = tcrs$epitope)
```

## Distance Distribution

Histogram and density of pairwise distances:

``` r
# Basic distribution plot (histogram + density + median line)
plot_distance_distribution(tcrs, organism = "mouse")

# Adjust number of bins
plot_distance_distribution(tcrs, organism = "mouse", bins = 50)
```

## CDR3 Sequence Logos

Display amino acid frequency at each position as a sequence logo:

``` r
cdr3_seqs <- tcrs$cdr3b

# Information content ("bits") mode
plot_cdr3_logo(cdr3_seqs, chain = "beta", method = "bits")

# Frequency ("prob") mode
plot_cdr3_logo(cdr3_seqs, chain = "beta", method = "prob")
```

The logo construction aligns sequences using BLOSUM62-optimal gap
placement, selects a center sequence, and builds a position weight
matrix (PWM) from the aligned set.

## Junction Composition Bars

Stacked bar chart showing V, N-insertion, D, and J nucleotide
composition at each CDR3 position:

``` r
# Requires nucleotide sequences
cdr3b_nucseqs <- c(
  "tgtgccagcagtgatttcgactggggaggggatgcagaaacgctgtatttt",
  "tgtgctagcagtccggacaggggtgaagtcttcttt",
  "tgtgctagcacagggggaggggctccgcttttt",
  "tgtgccagcagccaagatcctggggactatgaacagtacttc",
  "tgtgctagcagtcttggaggggaaaacacgctgtacttt"
)

plot_junction_bars(cdr3b_nucseqs, chain = "beta")
```

## Gene Usage

Horizontal bar chart of V-gene frequencies:

``` r
# V-alpha gene usage
plot_gene_usage(tcrs$va, chain = "alpha")

# V-beta gene usage
plot_gene_usage(tcrs$vb, chain = "beta")
```

## PCA Scatter Plot

2D scatter plot from kernel PCA or any dimensionality reduction:

``` r
pca <- compute_tcrdist_kernel_pca(tcrs, organism = "mouse",
                                   n_components = 5)

# Colored by epitope
plot_tcr_scatter(
  pca$embeddings[, 1], pca$embeddings[, 2],
  color = tcrs$epitope[pca$indices],
  xlab = "KPC1", ylab = "KPC2"
)

# Colored by a continuous variable
plot_tcr_scatter(
  pca$embeddings[, 1], pca$embeddings[, 2],
  color = pca$embeddings[, 3],
  xlab = "KPC1", ylab = "KPC2"
)
```

## Combining Plots with patchwork

Use the patchwork package to arrange multiple plots:

``` r
library(patchwork)

p1 <- plot_tcrdist_heatmap(tcrs, organism = "mouse")
p2 <- plot_distance_distribution(tcrs, organism = "mouse")
p3 <- plot_gene_usage(tcrs$va, chain = "alpha")
p4 <- plot_gene_usage(tcrs$vb, chain = "beta")

(p1 | p2) / (p3 | p4) +
  plot_annotation(title = "TCR Repertoire Overview")
```
