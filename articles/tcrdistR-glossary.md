# Glossary of Terms

This glossary defines computational, statistical, and biological terms
used throughout the tcrdistR documentation. It is organized by topic for
quick reference.

## TCR Biology

- **AIRR**: Adaptive Immune Receptor Repertoire. A community standard
  format for storing immune receptor sequencing data as tab-separated
  files.
- **CDR1, CDR2**: Complementarity-Determining Regions 1 and 2. Encoded
  entirely by the germline V-gene segment. They primarily contact the
  MHC molecule rather than the peptide antigen.
- **CDR2.5 (HV4)**: A tcrdist-specific term for the fourth hypervariable
  loop between CDR2 and CDR3. This region makes additional MHC contacts
  and is included in the TCRdist V-region distance.
- **CDR3**: Complementarity-Determining Region 3. The most variable part
  of the TCR, formed by V(D)J recombination. It is the primary
  determinant of antigen specificity — the amino acid sequence in this
  loop directly contacts the peptide antigen.
- **Clone / Clonotype**: A unique TCR sequence defined by its V-gene and
  CDR3 combination. A “clone” refers to a T cell or population sharing
  that sequence; a “clonotype” is the sequence identity itself.
- **Meta-clonotype**: A group of similar TCRs (within a distance
  threshold) found across multiple individuals. These represent
  convergent immune responses where different people independently
  generate similar TCRs against the same antigen. Meta-clonotypes can
  serve as biomarkers of antigen exposure.
- **Public clonotype**: A TCR sequence shared across multiple unrelated
  individuals, suggesting convergent selection driven by a common
  antigen.
- **Repertoire**: The full collection of TCR sequences in a biological
  sample.
- **V-region**: The portion of the TCR encoded by the V gene segment,
  encompassing the CDR1, CDR2, and CDR2.5 loops. These germline-encoded
  regions are less variable than CDR3 and primarily contact the MHC
  molecule.

## Distance and Sequence Analysis

- **BLOSUM62**: BLOcks SUbstitution Matrix, version 62. A widely used
  scoring matrix from protein evolution research that quantifies how
  often one amino acid substitutes for another in related proteins.
  Biochemically similar amino acids get high scores.
- **BSD4 matrix**: A BLOSUM62-derived substitution matrix rescaled
  specifically for TCR distance calculations (from Dash et al., 2017).
  This is the default scoring matrix used by TCRdist.
- **Gap penalty**: When aligning two CDR3 sequences of different
  lengths, a gap (insertion or deletion) must be introduced at one or
  more positions. The gap penalty (default: 12) is the distance cost
  added for each gap position, reflecting that insertions/deletions are
  rarer than amino acid substitutions.
- **Hamming distance**: The simplest sequence comparison: count the
  number of positions where two equal-length sequences have different
  characters.
- **Pairwise distance**: The distance computed between every possible
  pair in a set. For N TCRs, this produces an N x N symmetric matrix
  where entry (i, j) is the distance between TCR i and TCR j.
- **Substitution matrix**: A table scoring the similarity between all
  pairs of amino acids. High scores indicate biochemically similar amino
  acids that frequently interchange during evolution (e.g., leucine and
  isoleucine).
- **TCRdist**: The specific distance metric used by this package. It
  sums: (1) the CDR3 alignment distance (using BSD4 matrix with optimal
  gap placement), weighted 3x, plus (2) the V-region distance from CDR1,
  CDR2, and CDR2.5 loops. Lower distance = more similar TCRs. Distances
  below ~50 often indicate TCRs that recognize the same antigen.

## Statistics and Diversity

- **Clonality**: A measure of how dominated a repertoire is by a few
  clonotypes. Calculated as 1 - (Shannon entropy / log(number of
  clonotypes)). Ranges from 0 (all clonotypes equally abundant) to 1
  (single dominant clone).

- **Effective number of species**: The number of equally-abundant
  clonotypes that would give the same diversity index value. Converts an
  abstract diversity score into an intuitive count. For example, an
  effective number of 50 means the repertoire behaves as if it had 50
  equally-abundant clonotypes.

- **Gini index**: A measure of inequality in clonotype abundances,
  borrowed from economics. 0 = all clonotypes equally abundant (perfect
  equality); 1 = one clonotype has all the cells (maximum inequality).

- **Hill numbers**:

  A family of diversity indices parameterized by an order *q*. q=0 gives
  species richness (number of clonotypes), q=1 gives Shannon diversity,
  q=2 gives Simpson’s diversity. Higher orders give more weight to
  abundant clonotypes.

- **Jaccard index**: The fraction of clonotypes shared between two
  samples: \|intersection\| / \|union\|. Ignores abundance; only
  considers presence/absence. Ranges from 0 (no shared clonotypes) to 1
  (identical sets).

- **Morisita-Horn index**: An abundance-weighted overlap measure between
  two samples. Unlike Jaccard, it accounts for how many copies of each
  clonotype are present. Two samples sharing a single dominant clonotype
  will have a high Morisita-Horn value even if they differ in rare
  clonotypes.

- **Odds ratio**: In neighborhood tests, the ratio of odds of being the
  target category inside versus outside a TCR’s neighborhood. Values \>
  1 indicate enrichment; values \< 1 indicate depletion.

- **p-value (adjusted)**: The probability of observing a result this
  extreme by chance, corrected for the fact that many tests were
  performed simultaneously (one per TCR). tcrdistR uses the
  Benjamini-Hochberg correction, which controls the false discovery
  rate. A threshold of 0.05 means no more than 5% of significant results
  are expected to be false positives.

- **Shannon entropy**: An information-theoretic measure of diversity: H
  = -sum(p \* log(p)) where p is the proportion of each clonotype.
  Higher entropy = more diverse. For S equally-abundant clonotypes, H =
  log(S).

- **Simpson’s diversity**: The probability that two randomly chosen TCRs
  from the sample are different clonotypes. Higher values = more
  diverse. Ranges from 0 to 1.

- **Species richness**: The count of distinct clonotypes in a sample.
  The simplest diversity measure — it ignores abundances entirely.

## Dimensionality Reduction and Visualization

- **Dense matrix**: A standard matrix that stores every value. An N x N
  distance matrix uses approximately 8 \* N^2 bytes of memory (e.g.,
  ~800 MB for 10,000 TCRs).

- **Eigenvalues**: In PCA, eigenvalues quantify how much variance each
  principal component captures. Larger eigenvalues indicate more
  important components. Examining the dropoff in eigenvalues helps
  decide how many components are meaningful.

- **Kernel PCA**: A variant of PCA that works on a similarity (kernel)
  matrix rather than raw data coordinates. Used here to embed TCRdist
  distances into a low-dimensional space for visualization. “Linear
  kernel” applies a simple distance-to-similarity conversion; “Gaussian
  kernel” uses an exponential transformation that emphasizes local
  structure.

- **PCA (Principal Component Analysis)**: A method to reduce
  high-dimensional data to a few key axes (components) that capture the
  most variation. Points close together in PCA space have similar TCR
  sequences. The first few components typically explain the most
  biologically meaningful differences.

- **Sparse matrix**:

  A matrix that only stores non-zero (or below-threshold) entries. For
  distance matrices, entries above a threshold are not stored,
  dramatically reducing memory usage. Stored as a `dgCMatrix` object
  from the Matrix package.

- **UMAP (Uniform Manifold Approximation and Projection)**: A nonlinear
  dimensionality reduction method that often separates clusters better
  than PCA. Points close together in UMAP space share similar TCR
  sequences. Unlike PCA, UMAP preserves local neighborhood structure but
  does not preserve global distances.

## R Programming

- **S4 object**:

  An R class system (used by Bioconductor packages). tcrdistR’s `TCRrep`
  is an S4 object. Access its internal data using `@` instead of `$`.
  For example, `rep@clone_df` accesses the clonotype table stored inside
  a TCRrep object.

- **Faceting**:

  A ggplot2 concept: splitting a single plot into multiple panels based
  on a grouping variable (e.g., one panel per epitope).
  [`facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  creates a wrapped grid;
  [`facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html)
  creates a row-by-column grid.
