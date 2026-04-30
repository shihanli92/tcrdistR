# Advanced TCR Repertoire Analysis

## Setup

This vignette demonstrates the advanced analysis features in tcrdistR:
hierarchical clustering, neighborhood statistical tests, clumping
detection, meta-clonotype discovery, database matching, kernel PCA, and
more.

We use the DASH dataset — 1924 paired alpha-beta mouse TCRs across 7
viral epitopes from 78 subjects.

``` r
library(tcrdistR)
data(dash)

# Work with PA-specific TCRs for clustering demos (N = 324)
pa <- dash[dash$epitope == "PA", ]
nrow(pa)
#> [1] 324
```

## Hierarchical Clustering

Cluster TCRs by distance using hierarchical agglomerative clustering:

``` r
hc <- tcrdist_hclust(pa, organism = "mouse", method = "average")

# hc contains:
#   hc$hclust      - an hclust object (pass to plot(), cutree(), etc.)
#   hc$dist_matrix - the underlying distance matrix
#   hc$indices     - row indices used (relevant if subsampled)
names(hc)
#> [1] "hclust"      "dist_matrix" "indices"
```

Cut the dendrogram into a fixed number of clusters:

``` r
clusters <- cluster_tcrs(pa, organism = "mouse", k = 5)
table(clusters)
#> clusters
#>   1   2   3   4   5 
#>   9 265  45   2   3

# Or cut by height threshold
clusters_h <- cluster_tcrs(pa, organism = "mouse", h = 100)
table(clusters_h)
#> clusters_h
#>   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
#>   2   9   2   1  12   1  16   1  12   1   1   1   1  25   2  61   5   1   2   1 
#>  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40 
#>   1   1  20  15   4   2   2   7   2   1   6   4   1   1   1   1   1   1   4   3 
#>  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60 
#>   1   1   1   1   2   1   1   1   2   1   1   1   1   1   1   1   1   5   1   1 
#>  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80 
#>   1   3   1   1   1   1   1   1   1   1   1   2   1   1   1   1   1   2   1   1 
#>  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 
#>   1   1   1   1   1   1   2   1   2   1   1   1   1   1   1   1   1   1   2   1 
#> 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 
#>   1   1   2   1   1   1   1   1   1   1   1   1   1   1   1
```

## Neighborhood Statistical Tests

Test whether epitope specificity is enriched in TCR neighborhoods. This
computes TCRdist radius-based neighborhoods and runs a statistical test
for each TCR.

``` r
# Use a mixed-epitope subset for a meaningful test
mixed <- dash[dash$epitope %in% c("PA", "NP"), ]
mixed$is_PA <- ifelse(mixed$epitope == "PA", "PA", "other")

result <- neighborhood_test(
  mixed, organism = "mouse",
  variable = mixed$is_PA,
  radius = 50, test = "fisher"
)

head(result[order(result$p_adjusted), c("index", "n_neighbors", "odds_ratio",
                                         "p_value", "p_adjusted")])
#>     index n_neighbors odds_ratio     p_value   p_adjusted
#> 52     52          42        Inf 1.31689e-14 2.123907e-13
#> 72     72          42        Inf 1.31689e-14 2.123907e-13
#> 170   170          42        Inf 1.31689e-14 2.123907e-13
#> 171   171          42        Inf 1.31689e-14 2.123907e-13
#> 174   174          42        Inf 1.31689e-14 2.123907e-13
#> 175   175          42        Inf 1.31689e-14 2.123907e-13
```

For multi-level categorical variables, use the chi-squared test:

``` r
result_chi <- neighborhood_test(
  dash, organism = "mouse",
  variable = dash$epitope,
  radius = 50, test = "chisq"
)
```

## Meta-Clonotype Discovery

Identify quasi-public TCR motifs shared across multiple subjects. A
meta-clonotype is defined by a center TCR and a radius: all TCRs within
the radius from multiple subjects form the meta-clonotype.

``` r
meta <- find_meta_clonotypes(
  pa, organism = "mouse",
  radius = 48,
  min_nsubject = 2L,
  subject_col = "subject"
)

nrow(meta)
#> [1] 113
head(meta[, c("cdr3a", "cdr3b", "radius", "K_neighbors", "nsubject")])
#>              cdr3a          cdr3b radius K_neighbors nsubject
#> 1 CALGDRATGGNNKLTF   CASSPDRGEVFF     48           3        3
#> 2 CILRVGATGGNNKLTL   CASSLDRGEVFF     48          15       10
#> 3    CALVPSNTNKVVF   CASSLSGYEQYF     48          10        9
#> 4    CALGGGSNYQLIW    CASSLGGEVFF     48          29       12
#> 5  CVLSARAEGADRLTF CASSQAGDSYEQYF     48           4        4
#> 6 CAASGGTTASLGKLQF   CTCSADENTLYF     48           2        2
```

Summarize the composition of a discovered meta-clonotype:

``` r
if (nrow(meta) > 0) {
  indices <- as.integer(strsplit(meta$neighbor_indices[1], ",")[[1]])
  summary <- summarize_meta_clonotype(pa, meta$center_index[1], indices)
  summary$n_members
  summary$va_usage
  summary$vb_usage
}
#> TRBV29*01 
#>         3
```

## Clumping Detection

Test whether TCRs form statistically significant clusters beyond what is
expected by chance. This uses a background resampling model to compute
expected neighbor counts.

The
[`find_clumping()`](https://shihanli92.github.io/tcrdistR/reference/find_clumping.md)
function requires nucleotide sequences and J-gene annotations, which the
DASH dataset provides:

``` r
clump <- find_clumping(
  pa, organism = "mouse",
  radii = c(24L, 48L, 72L),
  num_random_samples = 50000L,
  verbose = TRUE
)

# How many TCRs are in significant clumps?
sum(clump$is_clumped)

# Detailed results sorted by significance
head(clump$results_df[, c("cdr3a", "cdr3b", "nbr_radius",
                            "pvalue_adj", "num_nbrs")])
```

## Database Matching

Match query TCRs against the built-in literature epitope database. This
requires background sampling and is available for human TCRs:

``` r
# Paired-chain matching with background-corrected p-values
results <- match_tcrs_to_db(human_tcrs, organism = "human")
significant <- results[results$pvalue_adj < 0.05, ]

# Single-chain exact CDR3 matching (no background model needed)
hits <- strict_single_chain_match_tcrs_to_db(human_tcrs, organism = "human")
hits$alpha_matches
hits$beta_matches
```

## Kernel PCA

Embed TCRs in lower-dimensional space via kernel PCA on the distance
matrix. The default kernel is linear (`1 - D/Dmax`); a Gaussian kernel
is also available.

``` r
pca <- compute_tcrdist_kernel_pca(
  pa, organism = "mouse",
  n_components = 10L,
  method = "eigen"
)

# Return structure
names(pca)
#> [1] "embeddings"   "eigenvalues"  "n_components"
dim(pca$embeddings)     # N x n_components
#> [1] 324  10
length(pca$eigenvalues)
#> [1] 10

# Top eigenvalues (variance captured)
head(pca$eigenvalues)
#> [1] 24.786698 15.525856 11.676165 10.940108  7.882001  5.825661
```

With Gaussian kernel:

``` r
pca_gauss <- compute_tcrdist_kernel_pca(
  pa, organism = "mouse",
  n_components = 10L,
  kernel = "gaussian",
  gaussian_kernel_sdev = 100
)
head(pca_gauss$eigenvalues)
#> [1] 35.031673 21.435992 16.936445 15.145947 11.507637  8.617358
```

## Fuzzy Diversity

Compute TCR-aware diversity that accounts for sequence similarity.
Unlike standard diversity which treats each unique clonotype as
completely distinct, fuzzy diversity merges clonotypes within a distance
threshold:

``` r
fd <- tcr_fuzzy_diversity(pa, organism = "mouse", threshold = 50)
fd$fuzzy_diversity     # diversity accounting for similar TCRs
#> [1] 0.977633
fd$standard_diversity  # standard Simpson's for comparison
#> [1] 1
```

## Distance-Based Joins

Fuzzy-join two TCR datasets by distance threshold. This finds all pairs
of TCRs (one from each dataset) within a specified radius:

``` r
left <- pa[1:20, ]
right <- pa[21:100, ]

# Inner join: only matched pairs within radius
joined <- tcrdist_join(left, right, organism = "mouse",
                       radius = 100, max_n = 3)
dim(joined)
#> [1] 24 25
head(joined[, c("cdr3a_x", "cdr3b_x", "cdr3a_y", "cdr3b_y", "tcrdist")])
#>            cdr3a_x       cdr3b_x          cdr3a_y      cdr3b_y tcrdist
#> 1 CALGDRATGGNNKLTF  CASSPDRGEVFF CALGGGATGGNNKLTF CASTPDRGEVFF      33
#> 2 CALGDRATGGNNKLTF  CASSPDRGEVFF  CALGDRRTGNYKYVF CASSLDRGEQYF      72
#> 3 CALGDRATGGNNKLTF  CASSPDRGEVFF    CALGAGSGNKLIF CASSWDRGEVFF      84
#> 4    CALVPSNTNKVVF CASSLGGENTLYF    CALAPSNTNKVVF CASSYGGGEQYF      72
#> 5    CALVPSNTNKVVF CASSLGGENTLYF    CALVPSNTNKVVF   CASSQAEVFF      81
#> 6   CALGKNYNQGKLIF   CASDRAGEQYF    CALGEGSNYKLTF  CASSLGGEQYF      93
```

## UMAP Visualization

Embed TCRs in 2D using UMAP computed directly from TCRdist KNN
neighbors:

``` r
umap <- compute_tcrdist_umap(pa, organism = "mouse", K = 15L, seed = 42)

# umap$embeddings is an N x 2 matrix
plot_tcr_scatter(
  umap$embeddings,
  color_by = pa$subject,
  axis_label_prefix = "UMAP",
  point_size = 2
)
```

## TCR Network Visualization

Build an igraph-based TCR similarity network. TCRs become nodes, edges
connect TCRs within a distance threshold. The threshold can be set
manually or auto-detected from the bimodal distance distribution:

``` r
# Auto-detect threshold from distance distribution
net <- compute_tcr_network(pa, organism = "mouse")

# Color by vertex attribute name (metadata stored in the graph)
plot_tcr_network(net, color_by = "epitope", vertex_size = 3)

# Manual threshold with minimum edge pruning
net2 <- compute_tcr_network(pa, organism = "mouse", threshold = 48,
                            min_edges = 2)
```

## CD8 Scoring

Score TCRs for CD8 phenotype using a logistic regression model trained
on known CD4/CD8 repertoires. This is only available for human TCRs:

``` r
scores <- make_cd8_score_table_column(human_tcrs)
# Returns a numeric vector of CD8 probability scores
```
