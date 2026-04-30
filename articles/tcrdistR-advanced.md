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

The neighborhood test asks a key biological question: **is antigen
specificity encoded in TCR sequence similarity?** For each TCR, it
defines a “neighborhood” (all TCRs within a distance radius), then tests
whether a variable of interest (e.g., epitope label) is enriched in the
neighborhood compared to the overall dataset.

A significant result (low adjusted p-value, high odds ratio) means that
knowing a TCR’s sequence tells you something about its specificity —
nearby TCRs tend to recognize the same antigen. The p-values are
adjusted for multiple testing using the Benjamini-Hochberg procedure
(controls false discovery rate).

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

Meta-clonotypes represent convergent immune responses: TCR motifs that
appear in multiple individuals, suggesting a shared solution to
recognizing the same antigen. These can serve as biomarkers of antigen
exposure, infection, or vaccine response.

A meta-clonotype is defined by a center TCR and a distance radius. All
TCRs within this radius from at least `min_nsubject` different subjects
form the meta-clonotype.

``` r
meta <- find_meta_clonotypes(
  pa, organism = "mouse",
  radius = 48,
  min_nsubject = 2L,
  subject_col = "subject"
)

nrow(meta)
#> [1] 324
head(meta[, c("cdr3a", "cdr3b", "radius", "K_neighbors", "nsubject")])
#>              cdr3a             cdr3b radius K_neighbors nsubject
#> 1    CAVSLDSNYQLIW CASSDFDWGGDAETLYF     48         324       15
#> 2 CALGDRATGGNNKLTF      CASSPDRGEVFF     48         324       15
#> 3   CALGSNTGYQNFYF       CASTGGGAPLF     48         324       15
#> 4    CALAPSNTNKVVF    CASSQDPGDYEQYF     48         324       15
#> 5    CALVPSNTNKVVF     CASSLGGENTLYF     48         324       15
#> 6   CALGKNYNQGKLIF       CASDRAGEQYF     48         324       15
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
#> 
#>   TRBV29*01   TRBV19*01   TRBV19*03 TRBV13-3*01    TRBV2*01   TRBV17*01 
#>         209          32          19          12          12          11 
#> TRBV13-1*01   TRBV14*01    TRBV5*01    TRBV1*01 TRBV13-2*01   TRBV31*01 
#>           5           4           4           3           3           3 
#>    TRBV4*01   TRBV20*01   TRBV15*01   TRBV24*02 
#>           3           2           1           1
```

## Clumping Detection

Clumping asks whether TCR neighborhoods are larger than expected by
chance. Even in a random repertoire, some TCRs will be close to each
other simply because they use the same V-gene or have similar CDR3
lengths. The clumping test corrects for this by building a background
model through chain resampling (shuffling alpha and beta chains
independently to break biologically meaningful pairing) and comparing
observed neighbor counts against this null distribution using a Poisson
model (Dash et al., 2017).

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

### Human TCR example with the flu dataset

tcrdistR also ships with a **flu** dataset containing 2271 human TCRs
from VDJdb. This dataset works identically with `organism = "human"`:

``` r
data(flu)
flu_sub <- flu[1:200, ]

pca_flu <- compute_tcrdist_kernel_pca(
  flu_sub, organism = "human", n_components = 10L
)
dim(pca_flu$embeddings)
#> [1] 200  10
```

## Fuzzy Diversity

Standard diversity treats each unique CDR3 sequence as completely
distinct, even if two sequences differ by a single amino acid. Fuzzy
diversity accounts for sequence similarity: clonotypes within a TCRdist
threshold are considered “the same” for diversity calculations. This
gives a more biologically meaningful measure when a repertoire contains
many closely related variants (e.g., from somatic hypermutation or
convergent recombination):

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

## Session Info

``` r
sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] tcrdistR_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6         knitr_1.51        rlang_1.2.0       xfun_0.57        
#>  [5] textshaping_1.0.5 jsonlite_2.0.0    htmltools_0.5.9   ragg_1.5.2       
#>  [9] sass_0.4.10       rmarkdown_2.31    grid_4.6.0        evaluate_1.0.5   
#> [13] jquerylib_0.1.4   fastmap_1.2.0     yaml_2.3.12       lifecycle_1.0.5  
#> [17] compiler_4.6.0    RSpectra_0.16-2   fs_2.1.0          Rcpp_1.1.1-1.1   
#> [21] systemfonts_1.3.2 lattice_0.22-9    digest_0.6.39     R6_2.6.1         
#> [25] bslib_0.10.0      Matrix_1.7-5      tools_4.6.0       pkgdown_2.2.0    
#> [29] cachem_1.1.0      desc_1.4.3
```
