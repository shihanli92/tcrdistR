# tcrdistR Development Progress

## Project Overview

Standalone R package for TCRdist — TCR distance calculations and
repertoire analysis. Targeting feature parity with Python `tcrdist3`.
C++-first architecture via Rcpp.

**Source reference:** `rconga` package (existing TCRdist R+C++ code to
extract and extend) **Python reference:** `tcrdist3` (kmayerb/tcrdist3)
for full feature set

------------------------------------------------------------------------

## Architecture: C++-First

All computation in C++ via Rcpp. R layer is thin wrappers (input
validation, S4 class glue, I/O) and plotting only.

------------------------------------------------------------------------

## Current State

**992 tests, R CMD check: 0 errors, 0 warnings, 6 NOTEs** (C++14,
extdata size, etc.)

### Completed Phases

#### Phase 1: Foundation — C++ Core + Minimal R Shell ✅

- `src/tcrdist_core.h` — shared AA index, CDR3Data, cdr3_dist_fast,
  VDistLookup, budgeted overload
- `src/blosum.cpp` — BLOSUM62/BSD4 in C++ (constexpr)
- `src/tcrdist_distances.cpp` — single-pair CDR3 + V-region + paired
- `src/v_region_distances.cpp` — V-gene pairwise distances in C++
- `src/tcrdist_matrix.cpp` — dense N×N
- `R/all_genes.R`, `R/constants.R`, `R/wrappers_distance.R`,
  `R/tcrdistR-package.R`
- `inst/extdata/combo_xcr_2023-12-30.tsv` — gene database

#### Phase 2: S4 Classes + Sparse/Rectangular + KNN ✅

- `R/classes.R` — S4 TCRrep class + `R/constructors.R` +
  `R/show-methods.R`
- `src/tcrdist_sparse.cpp` + `R/wrappers_sparse.R`
- `src/tcrdist_rect.cpp` + `R/wrappers_rectangular.R`
- `src/tcrdist_knn.cpp` + `src/tcrdist_radius.cpp` +
  `R/wrappers_neighbors.R`
- `src/find_neighbors.cpp` — KNN from precomputed matrices + PCA
  embeddings

#### Phase 3: I/O Parsers ✅

- `R/io_airr.R`, `R/io_adaptive.R`, `R/io_10x.R`, `R/io_generic.R`
- `R/io_utils.R` — shared gene name normalization
- `R/translation.R` — nucleotide-to-protein translation, genetic code

#### Phase 4: Clumping + Background Models ✅

- `src/tcrdist_background.cpp` — background frequency computation
- `src/tcr_clumping.cpp` — Poisson test loop
- `src/tcr_sampler.cpp` — junction parsing, allele optimization, chain
  resampling
- `R/tcr_clumping.R` — find_clumping(), setup_tcr_groups(),
  single-linkage clustering
- `R/tcr_sampler.R` — junction analysis, allele finding, resampling
  wrappers

#### Phase 5: DB Matching, Kernel PCA, CD8 Scoring ✅

- `R/tcr_db_matching.R` — find_significant_tcrdist_matches(),
  match_tcrs_to_db(), strict_single_chain_match_tcrs_to_db()
- `R/kernel_pca.R` — compute_tcrdist_kernel_pca(), scipy-validated
  (eigenvalues to 1e-10, embeddings to 1e-8)
- `R/cd8_scoring.R` — make_cd8_score_table_column(), logistic regression
  with cached model loading
- `inst/extdata/` — TCR literature databases (human paired, human/mouse
  single-chain), CD8 model weights
- `tests/generate_kernel_pca_fixture.py` — Python script to regenerate
  scipy reference
- `tests/testthat/fixtures/dash.csv`, `kernel_pca_ref.json` —
  cross-validation fixtures

------------------------------------------------------------------------

## Package Structure

    tcrdistR/
    ├── R/                            # 36 R files (~9,900 lines)
    │   ├── tcrdistR-package.R        # Package doc, .onLoad, .tcrdistR_env
    │   ├── classes.R                 # S4 class: TCRrep
    │   ├── constructors.R            # TCRrep() constructor (with clone deduplication)
    │   ├── show-methods.R            # print/show/summary
    │   ├── constants.R               # AMINO_ACIDS, weights, gap penalties, BLOSUM62
    │   ├── all_genes.R               # Gene database loading, V-distance matrices
    │   ├── genetic_code.R            # Genetic code constants
    │   ├── translation.R             # Nucleotide translation
    │   ├── wrappers_distance.R       # tcrdist_matrix(), hamming_distance/matrix(), .resolve_components()
    │   ├── wrappers_neighbors.R      # tcrdist_knn(), tcrdist_radius_neighbors(), knn_from_matrix/pca
    │   ├── wrappers_sparse.R         # tcrdist_sparse()
    │   ├── wrappers_rectangular.R    # tcrdist_rect()
    │   ├── diversity.R               # tcr_diversity(), tcr_fuzzy_diversity(), tcr_richness(), tcr_clonality()
    │   ├── hierarchical.R            # tcrdist_hclust(), cluster_tcrs(), neighborhood_test()
    │   ├── meta_clonotypes.R         # find_meta_clonotypes(), summarize_meta_clonotype()
    │   ├── tcr_join.R                # tcrdist_join() distance-based fuzzy joins
    │   ├── tcr_clumping.R            # find_clumping(), setup_tcr_groups()
    │   ├── tcr_sampler.R             # Junction parsing, allele optimization, resampling
    │   ├── tcr_db_matching.R         # DB matching with background-corrected p-values
    │   ├── kernel_pca.R              # Kernel PCA (scipy-validated)
    │   ├── cd8_scoring.R             # CD8 logistic regression scoring
    │   ├── umap.R                    # compute_tcrdist_umap() via rconga KNN
    │   ├── network.R                 # compute_tcr_network(), plot_tcr_network() (igraph)
    │   ├── as_tcr_df.R               # as_tcr_df() data standardization
    │   ├── data.R                    # Data object documentation
    │   ├── io_airr.R                 # AIRR format reader
    │   ├── io_adaptive.R             # Adaptive ImmunoSeq reader
    │   ├── io_10x.R                  # 10x Genomics reader
    │   ├── io_generic.R              # Generic CSV/TSV reader
    │   ├── io_utils.R                # Shared I/O utilities
    │   ├── plot_utils.R              # Palette helpers, ggplot2 check
    │   ├── plot_logo.R               # CDR3 logos, junction bars, gene usage
    │   ├── plot_distances.R          # Heatmap, dendrogram, distribution
    │   ├── plot_scatter.R            # Generic 2D scatter (PCA/UMAP)
    │   └── RcppExports.R             # Auto-generated
    ├── src/                          # 13 C++ files (~3,100 lines)
    │   ├── tcrdist_core.h            # Shared: CDR3Data, cdr3_dist_fast, VDistLookup, BSD4
    │   ├── blosum.cpp                # BLOSUM62/BSD4 (constexpr)
    │   ├── tcrdist_distances.cpp     # Single-pair distances
    │   ├── tcrdist_matrix.cpp        # Dense N×N pairwise matrix
    │   ├── tcrdist_sparse.cpp        # Sparse matrix + early termination
    │   ├── tcrdist_rect.cpp          # Rectangular query-vs-reference
    │   ├── tcrdist_knn.cpp           # K-nearest neighbors
    │   ├── tcrdist_radius.cpp        # Radius-based neighbors
    │   ├── v_region_distances.cpp    # V-gene pairwise distances
    │   ├── find_neighbors.cpp        # KNN from precomputed/PCA matrices
    │   ├── tcrdist_background.cpp    # Background distributions
    │   ├── tcr_clumping.cpp          # Poisson test loop
    │   ├── tcr_sampler.cpp           # Junction parsing, resampling
    │   ├── tcrdist_hamming.cpp       # Hamming distance (single-pair + matrix)
    │   └── RcppExports.cpp           # Auto-generated
    ├── inst/extdata/
    │   ├── combo_xcr_2023-12-30.tsv            # Gene database (1.2 MB)
    │   ├── new_paired_tcr_db_for_matching_nr.tsv  # Paired TCR DB, human (507 KB)
    │   ├── human_tcr_db_for_matching.tsv       # Single-chain DB, human (6.4 MB)
    │   ├── mouse_tcr_db_for_matching.tsv       # Single-chain DB, mouse (882 KB)
    │   ├── cd8_logreg_params_A.txt             # CD8 model, alpha chain (10 KB)
    │   └── cd8_logreg_params_B.txt             # CD8 model, beta chain (9.4 KB)
    ├── tests/testthat/               # 24 test files, 941 tests
    │   ├── test-blosum.R
    │   ├── test-cdr3-distance.R
    │   ├── test-tcrdist.R
    │   ├── test-gene-database.R
    │   ├── test-classes.R
    │   ├── test-sparse.R
    │   ├── test-rectangular.R
    │   ├── test-neighbors.R
    │   ├── test-io.R
    │   ├── test-clumping.R
    │   ├── test-db-matching.R
    │   ├── test-kernel-pca.R         # Cross-validated against scipy/sklearn
    │   ├── test-cd8-scoring.R
    │   ├── test-plotting.R           # Visualization tests (85 tests)
    │   ├── test-hamming.R            # Hamming distance tests
    │   ├── test-diversity.R          # Diversity metric tests
    │   ├── test-hierarchical.R       # Clustering + neighborhood tests
    │   ├── test-join.R               # Distance-based join tests
    │   ├── test-meta-clonotypes.R    # Meta-clonotype detection tests
    │   ├── test-components.R         # Per-component distance tests (54 tests)
    │   ├── test-network.R            # Network visualization tests
    │   ├── test-umap.R               # UMAP tests
    │   ├── test-dedup.R              # Clone deduplication tests
    │   └── fixtures/                 # dash.csv, kernel_pca_ref.json, I/O fixtures
    └── tests/generate_kernel_pca_fixture.py  # Regenerates scipy reference

------------------------------------------------------------------------

## Exported Functions (~70)

**Distance computation:** `tcrdist_matrix`, `tcrdist_rect`,
`tcrdist_sparse`, `weighted_cdr3_distance`, `bsd4_matrix`,
`hamming_distance`, `hamming_matrix` **Neighbor search:** `tcrdist_knn`,
`tcrdist_radius_neighbors`, `knn_from_matrix`, `knn_from_pca`
**Clumping:** `find_clumping`, `setup_tcr_groups` **DB matching:**
`find_significant_tcrdist_matches`, `match_tcrs_to_db`,
`strict_single_chain_match_tcrs_to_db` **Kernel PCA:**
`compute_tcrdist_kernel_pca` **UMAP:** `compute_tcrdist_umap` **CD8
scoring:** `make_cd8_score_table_column` **Diversity:** `tcr_diversity`,
`tcr_fuzzy_diversity`, `tcr_richness`, `tcr_clonality` **Clustering:**
`tcrdist_hclust`, `cluster_tcrs`, `neighborhood_test`
**Meta-clonotypes:** `find_meta_clonotypes`, `summarize_meta_clonotype`
**Joins:** `tcrdist_join` **Network:** `compute_tcr_network`,
`plot_tcr_network` **Visualization:** `plot_cdr3_logo`,
`plot_junction_bars`, `plot_gene_usage`, `plot_tcrdist_heatmap`,
`plot_tcrdist_dendrogram`, `plot_distance_distribution`,
`plot_tcr_scatter` **I/O:** `read_airr`, `read_adaptive`, `read_10x`,
`read_tcr_table`, `as_tcr_df` **Utilities:** `load_gene_database`,
`get_translation`, `reverse_complement`, `trim_allele_to_gene`,
`AMINO_ACIDS`, `TCRrep`

#### Phase 6: Visualization ✅

- `R/plot_utils.R` — `.tcrdistR_palette()`, `.check_ggplot2()` shared
  helpers
- `R/plot_logo.R` — `.align_cdr3_regions()`, `.build_cdr3_pwm()`,
  [`plot_cdr3_logo()`](https://shihanli92.github.io/tcrdistR/reference/plot_cdr3_logo.md),
  [`plot_junction_bars()`](https://shihanli92.github.io/tcrdistR/reference/plot_junction_bars.md),
  [`plot_gene_usage()`](https://shihanli92.github.io/tcrdistR/reference/plot_gene_usage.md)
- `R/plot_distances.R` — `.hclust_to_segments()`,
  [`plot_tcrdist_heatmap()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_heatmap.md),
  [`plot_tcrdist_dendrogram()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcrdist_dendrogram.md),
  [`plot_distance_distribution()`](https://shihanli92.github.io/tcrdistR/reference/plot_distance_distribution.md)
- `R/plot_scatter.R` —
  [`plot_tcr_scatter()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_scatter.md)
  (kernel PCA / UMAP 2D scatter)
- `R/constants.R` — added BLOSUM62 substitution matrix for CDR3
  alignment scoring
- Dependencies: ggplot2, ggseqlogo, patchwork, ggrepel (all Suggests)

#### Phase 7: Additional Features ✅

- `src/tcrdist_hamming.cpp` — Hamming distance (single-pair + N×N
  matrix) in C++
- `R/wrappers_distance.R` —
  [`hamming_distance()`](https://shihanli92.github.io/tcrdistR/reference/hamming_distance.md),
  [`hamming_matrix()`](https://shihanli92.github.io/tcrdistR/reference/hamming_matrix.md)
  R wrappers
- `R/diversity.R` —
  [`tcr_diversity()`](https://shihanli92.github.io/tcrdistR/reference/tcr_diversity.md)
  (generalized Simpson’s with CI),
  [`tcr_fuzzy_diversity()`](https://shihanli92.github.io/tcrdistR/reference/tcr_fuzzy_diversity.md)
  (TCR-aware diversity),
  [`tcr_richness()`](https://shihanli92.github.io/tcrdistR/reference/tcr_richness.md),
  [`tcr_clonality()`](https://shihanli92.github.io/tcrdistR/reference/tcr_clonality.md)
- `R/tcr_join.R` —
  [`tcrdist_join()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_join.md)
  distance-based fuzzy join (inner/left, max_n, suffix handling)
- `R/hierarchical.R` —
  [`tcrdist_hclust()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_hclust.md),
  [`cluster_tcrs()`](https://shihanli92.github.io/tcrdistR/reference/cluster_tcrs.md),
  [`neighborhood_test()`](https://shihanli92.github.io/tcrdistR/reference/neighborhood_test.md)
  (Fisher/chi-sq per TCR)
- `R/meta_clonotypes.R` —
  [`find_meta_clonotypes()`](https://shihanli92.github.io/tcrdistR/reference/find_meta_clonotypes.md)
  (ECDF-based radius, subject counting, subset de-dup),
  [`summarize_meta_clonotype()`](https://shihanli92.github.io/tcrdistR/reference/summarize_meta_clonotype.md)

#### Phase 8: Polish ✅

- `README.Rmd` / `README.md` — GitHub landing page with installation,
  quick start, feature overview
- `vignettes/tcrdistR-getting-started.Rmd` — Tutorial: data loading,
  distances, neighbors, diversity
- `vignettes/tcrdistR-advanced.Rmd` — Clustering, clumping,
  meta-clonotypes, DB matching, kernel PCA, joins
- `vignettes/tcrdistR-visualization.Rmd` — All 7 plot functions with
  examples
- `_pkgdown.yml` — Site config with 15 reference groups, 3 articles,
  navbar
- `@seealso` cross-references added to ~40 exported functions across 22
  R files
- `@examples` added to 5 previously-missing functions
- Package-level documentation (`tcrdistR-package.R`) rewritten with full
  narrative overview
- `benchmarks/` — R + Python benchmark scripts for tcrdistR vs tcrdist3
  comparison

#### Phase 9: Post-Release Enhancements ✅

- `R/constructors.R` — Clone deduplication in TCRrep constructor
  (matching tcrdist3 behavior)
- `vignettes/tcrdistR-tcrrep-workflow.Rmd` — End-to-end TCRrep workflow
  vignette
- `R/umap.R` —
  [`compute_tcrdist_umap()`](https://shihanli92.github.io/tcrdistR/reference/compute_tcrdist_umap.md)
  using rconga KNN pipeline
- `R/as_tcr_df.R` —
  [`as_tcr_df()`](https://shihanli92.github.io/tcrdistR/reference/as_tcr_df.md)
  for standardizing TCR data.frames from common tools (10x, Adaptive,
  AIRR, etc.)
- `R/network.R` —
  [`compute_tcr_network()`](https://shihanli92.github.io/tcrdistR/reference/compute_tcr_network.md),
  [`plot_tcr_network()`](https://shihanli92.github.io/tcrdistR/reference/plot_tcr_network.md)
  (igraph-based TCR networks with auto-threshold detection via bimodal
  KDE, node jittering, min_edges pruning)
- Per-component TCRdist distances — `components` parameter added to all
  5 distance wrappers (`tcrdist_matrix`, `tcrdist_sparse`,
  `tcrdist_rect`, `tcrdist_knn`, `tcrdist_radius_neighbors`); presets:
  `"all"`, `"cdr3"`, `"v_region"`, `"alpha"`, `"beta"`, plus custom
  combos. Per-chain CDR3 weight/gap params in all C++ entry points.
  Cross-validated against tcrdist3 on human (influenza) and mouse (DASH)
  datasets.
- `tests/testthat/test-components.R` — 54 component selection tests
  (additivity, consistency across wrappers)
- `tests/testthat/test-network.R` — 364-line network visualization test
  suite
- `tests/testthat/test-umap.R` — 190-line UMAP test suite
- `tests/testthat/test-dedup.R` — 200-line clone deduplication tests

#### Phase 10: Single-Chain Mode ✅

- [`.fill_missing_chain()`](https://shihanli92.github.io/tcrdistR/reference/dot-fill_missing_chain.md)
  helper detects alpha-only or beta-only input, fills dummy chain values
- All 5 distance wrappers (`tcrdist_matrix`, `tcrdist_sparse`,
  `tcrdist_rect`, `tcrdist_knn`, `tcrdist_radius_neighbors`) auto-detect
  single-chain input and set `components` accordingly
- [`tcrdist_rect()`](https://shihanli92.github.io/tcrdistR/reference/tcrdist_rect.md)
  validates that query and ref have matching chain types
- TCRrep constructor passes `components` based on `chains` slot when
  `compute_distances = TRUE`
- `tests/testthat/test-single-chain.R` — 50 tests covering all wrappers,
  TCRrep, error cases, component overrides, mouse organism

#### Phase 11: Documentation Polish ✅

- `vignettes/tcrdistR-getting-started.Rmd` — added sections:
  Per-Component Distances, Single-Chain Mode, Standardizing TCR Data
  (`as_tcr_df`)
- `vignettes/tcrdistR-advanced.Rmd` — added sections: UMAP
  Visualization, TCR Network Visualization
- `vignettes/tcrdistR-visualization.Rmd` — added section: TCR Network
  Plot
- `vignettes/tcrdistR-tcrrep-workflow.Rmd` — added section: Single-Chain
  TCRrep
- `README.Rmd` — added per-component and single-chain examples to quick
  start, updated description

------------------------------------------------------------------------

## Verification

1.  C++ outputs match known reference values from Python tcrdist3
2.  Kernel PCA matches scipy.linalg.eigh to 1e-10 (eigenvalues) / 1e-8
    (embeddings)
3.  Sparse/dense equivalence verified
4.  Rectangular correctness: rect(A,A) == matrix(A)
5.  R CMD check clean at every phase
6.  Cross-validation against rconga on same inputs
