# tcrdistR Development Progress

## Project Overview

Standalone R package for TCRdist — TCR distance calculations and repertoire analysis.
Targeting feature parity with Python `tcrdist3`. C++-first architecture via Rcpp.

**Source reference:** `rconga` package (existing TCRdist R+C++ code to extract and extend)
**Python reference:** `tcrdist3` (kmayerb/tcrdist3) for full feature set

---

## Architecture: C++-First

All computation in C++ via Rcpp. R layer is thin wrappers (input validation, S4 class glue, I/O) and plotting only.

---

## Current State

**583 tests, R CMD check: 0 errors, 0 warnings, 2 NOTEs** (C++14, extdata size)

### Completed Phases

#### Phase 1: Foundation — C++ Core + Minimal R Shell ✅
- `src/tcrdist_core.h` — shared AA index, CDR3Data, cdr3_dist_fast, VDistLookup, budgeted overload
- `src/blosum.cpp` — BLOSUM62/BSD4 in C++ (constexpr)
- `src/tcrdist_distances.cpp` — single-pair CDR3 + V-region + paired
- `src/v_region_distances.cpp` — V-gene pairwise distances in C++
- `src/tcrdist_matrix.cpp` — dense N×N
- `R/all_genes.R`, `R/constants.R`, `R/wrappers_distance.R`, `R/tcrdistR-package.R`
- `inst/extdata/combo_xcr_2023-12-30.tsv` — gene database

#### Phase 2: S4 Classes + Sparse/Rectangular + KNN ✅
- `R/classes.R` — S4 TCRrep class + `R/constructors.R` + `R/show-methods.R`
- `src/tcrdist_sparse.cpp` + `R/wrappers_sparse.R`
- `src/tcrdist_rect.cpp` + `R/wrappers_rectangular.R`
- `src/tcrdist_knn.cpp` + `src/tcrdist_radius.cpp` + `R/wrappers_neighbors.R`
- `src/find_neighbors.cpp` — KNN from precomputed matrices + PCA embeddings

#### Phase 3: I/O Parsers ✅
- `R/io_airr.R`, `R/io_adaptive.R`, `R/io_10x.R`, `R/io_generic.R`
- `R/io_utils.R` — shared gene name normalization
- `R/translation.R` — nucleotide-to-protein translation, genetic code

#### Phase 4: Clumping + Background Models ✅
- `src/tcrdist_background.cpp` — background frequency computation
- `src/tcr_clumping.cpp` — Poisson test loop
- `src/tcr_sampler.cpp` — junction parsing, allele optimization, chain resampling
- `R/tcr_clumping.R` — find_clumping(), setup_tcr_groups(), single-linkage clustering
- `R/tcr_sampler.R` — junction analysis, allele finding, resampling wrappers

#### Phase 5: DB Matching, Kernel PCA, CD8 Scoring ✅
- `R/tcr_db_matching.R` — find_significant_tcrdist_matches(), match_tcrs_to_db(), strict_single_chain_match_tcrs_to_db()
- `R/kernel_pca.R` — compute_tcrdist_kernel_pca(), scipy-validated (eigenvalues to 1e-10, embeddings to 1e-8)
- `R/cd8_scoring.R` — make_cd8_score_table_column(), logistic regression with cached model loading
- `inst/extdata/` — TCR literature databases (human paired, human/mouse single-chain), CD8 model weights
- `tests/generate_kernel_pca_fixture.py` — Python script to regenerate scipy reference
- `tests/testthat/fixtures/dash.csv`, `kernel_pca_ref.json` — cross-validation fixtures

---

## Package Structure

```
tcrdistR/
├── R/                            # 21 R files (~7,400 lines)
│   ├── tcrdistR-package.R        # Package doc, .onLoad, .tcrdistR_env
│   ├── classes.R                 # S4 class: TCRrep
│   ├── constructors.R            # TCRrep() constructor
│   ├── show-methods.R            # print/show/summary
│   ├── constants.R               # AMINO_ACIDS, weights, gap penalties
│   ├── all_genes.R               # Gene database loading, V-distance matrices
│   ├── genetic_code.R            # Genetic code constants
│   ├── translation.R             # Nucleotide translation
│   ├── wrappers_distance.R       # tcrdist_matrix(), weighted_cdr3_distance()
│   ├── wrappers_neighbors.R      # tcrdist_knn(), tcrdist_radius_neighbors(), knn_from_matrix/pca
│   ├── wrappers_sparse.R         # tcrdist_sparse()
│   ├── wrappers_rectangular.R    # tcrdist_rect()
│   ├── tcr_clumping.R            # find_clumping(), setup_tcr_groups()
│   ├── tcr_sampler.R             # Junction parsing, allele optimization, resampling
│   ├── tcr_db_matching.R         # DB matching with background-corrected p-values
│   ├── kernel_pca.R              # Kernel PCA (scipy-validated)
│   ├── cd8_scoring.R             # CD8 logistic regression scoring
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
├── src/                          # 12 C++ files (~3,000 lines)
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
│   └── RcppExports.cpp           # Auto-generated
├── inst/extdata/
│   ├── combo_xcr_2023-12-30.tsv            # Gene database (1.2 MB)
│   ├── new_paired_tcr_db_for_matching_nr.tsv  # Paired TCR DB, human (507 KB)
│   ├── human_tcr_db_for_matching.tsv       # Single-chain DB, human (6.4 MB)
│   ├── mouse_tcr_db_for_matching.tsv       # Single-chain DB, mouse (882 KB)
│   ├── cd8_logreg_params_A.txt             # CD8 model, alpha chain (10 KB)
│   └── cd8_logreg_params_B.txt             # CD8 model, beta chain (9.4 KB)
├── tests/testthat/               # 13 test files, 583 tests
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
│   └── fixtures/                 # dash.csv, kernel_pca_ref.json, I/O fixtures
└── tests/generate_kernel_pca_fixture.py  # Regenerates scipy reference
```

---

## Exported Functions (51)

**Distance computation:** `tcrdist_matrix`, `tcrdist_rect`, `tcrdist_sparse`, `weighted_cdr3_distance`, `bsd4_matrix`
**Neighbor search:** `tcrdist_knn`, `tcrdist_radius_neighbors`, `knn_from_matrix`, `knn_from_pca`
**Clumping:** `find_clumping`, `setup_tcr_groups`
**DB matching:** `find_significant_tcrdist_matches`, `match_tcrs_to_db`, `strict_single_chain_match_tcrs_to_db`
**Kernel PCA:** `compute_tcrdist_kernel_pca`
**CD8 scoring:** `make_cd8_score_table_column`
**Visualization:** `plot_cdr3_logo`, `plot_junction_bars`, `plot_gene_usage`, `plot_tcrdist_heatmap`, `plot_tcrdist_dendrogram`, `plot_distance_distribution`, `plot_tcr_scatter`
**I/O:** `read_airr`, `read_adaptive`, `read_10x`, `read_tcr_table`
**Utilities:** `load_gene_database`, `get_translation`, `reverse_complement`, `trim_allele_to_gene`, `AMINO_ACIDS`, `TCRrep`
**Rcpp (advanced):** 14 `rcpp_*` functions for direct C++ access

#### Phase 6: Visualization ✅
- `R/plot_utils.R` — `.tcrdistR_palette()`, `.check_ggplot2()` shared helpers
- `R/plot_logo.R` — `.align_cdr3_regions()`, `.build_cdr3_pwm()`, `plot_cdr3_logo()`, `plot_junction_bars()`, `plot_gene_usage()`
- `R/plot_distances.R` — `.hclust_to_segments()`, `plot_tcrdist_heatmap()`, `plot_tcrdist_dendrogram()`, `plot_distance_distribution()`
- `R/plot_scatter.R` — `plot_tcr_scatter()` (kernel PCA / UMAP 2D scatter)
- `R/constants.R` — added BLOSUM62 substitution matrix for CDR3 alignment scoring
- Dependencies: ggplot2, ggseqlogo, patchwork, ggrepel (all Suggests)

---

## Remaining Work

### Phase 7: Additional Features
- Meta-clonotype detection
- Additional distance metrics (Hamming, Levenshtein, Needleman-Wunsch)
- Hierarchical clustering
- Gene usage statistics
- TCR join (query matching)

### Phase 8: Polish
- Vignettes
- pkgdown site
- Benchmarks vs Python tcrdist3

---

## Verification

1. C++ outputs match known reference values from Python tcrdist3
2. Kernel PCA matches scipy.linalg.eigh to 1e-10 (eigenvalues) / 1e-8 (embeddings)
3. Sparse/dense equivalence verified
4. Rectangular correctness: rect(A,A) == matrix(A)
5. R CMD check clean at every phase
6. Cross-validation against rconga on same inputs
