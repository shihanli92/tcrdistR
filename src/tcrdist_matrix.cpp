// tcrdist_matrix.cpp — Full N x N pairwise TCRdist matrix.
//
// Exports:
//   rcpp_tcrdist_matrix() — N x N paired-chain distance matrix.
//
// Algorithm matches rconga/src/tcrdist_rcpp.cpp lines 198-370 exactly, but
// uses BSD4_FLAT (constexpr, no bsd4 matrix argument) and the shared
// CDR3Data / cdr3_dist_fast() infrastructure from tcrdist_core.h.

#include "tcrdist_core.h"

using namespace Rcpp;

//' Pairwise TCRdist matrix (C++ implementation)
//'
//' Computes the full N x N symmetric matrix of paired-chain TCRdist distances.
//' For each pair \code{(i, j)}:
//' \deqn{
//'   d(i,j) = v\_dist\_a[va_i, va_j] + cdr3\_dist(cdr3a_i, cdr3a_j)
//'           + v\_dist\_b[vb_i, vb_j] + cdr3\_dist(cdr3b_i, cdr3b_j)
//' }
//'
//' Only the upper triangle is computed; it is mirrored to the lower triangle.
//' The diagonal is zero.  \code{Rcpp::checkUserInterrupt()} is called every
//' 100 outer-loop rows to allow the user to interrupt long computations.
//'
//' BSD4 values are compiled in as \code{constexpr}; no \code{bsd4} matrix
//' argument is needed (unlike the rconga equivalent).
//'
//' @param va_genes   Character vector of length N. Alpha V-gene allele names.
//' @param cdr3a_seqs Character vector of length N. Alpha CDR3 sequences.
//' @param vb_genes   Character vector of length N. Beta V-gene allele names.
//' @param cdr3b_seqs Character vector of length N. Beta CDR3 sequences.
//' @param v_dist_a   Named square \code{NumericMatrix}. Pre-computed pairwise
//'   V-region distances for alpha genes (row/column names are V-gene allele
//'   strings, e.g. \code{"TRAV1-1*01"}).
//' @param v_dist_b   Named square \code{NumericMatrix}. Same structure for
//'   beta genes.
//' @param weight_cdr3_a Integer. Alpha CDR3 alignment distance multiplier.
//'   Default 3.
//' @param gap_penalty_cdr3_a Integer. Alpha per-residue length-difference
//'   penalty. Default 12.
//' @param weight_cdr3_b Integer. Beta CDR3 alignment distance multiplier.
//'   Default 3.
//' @param gap_penalty_cdr3_b Integer. Beta per-residue length-difference
//'   penalty. Default 12.
//' @return An N x N symmetric \code{NumericMatrix} of TCRdist distances.
//' @examples
//' \dontrun{
//'   mat <- rcpp_tcrdist_matrix(
//'     va_genes = c("TRAV1-1*01", "TRAV1-2*01"),
//'     cdr3a_seqs = c("CAVSANSGTYF", "CAVSANSGTYF"),
//'     vb_genes = c("TRBV20-1*01", "TRBV20-1*01"),
//'     cdr3b_seqs = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
//'     v_dist_a = v_alpha_mat,
//'     v_dist_b = v_beta_mat
//'   )
//' }
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rcpp_tcrdist_matrix(
    const CharacterVector& va_genes,
    const CharacterVector& cdr3a_seqs,
    const CharacterVector& vb_genes,
    const CharacterVector& cdr3b_seqs,
    const NumericMatrix&   v_dist_a,
    const NumericMatrix&   v_dist_b,
    int weight_cdr3_a       = 3,
    int gap_penalty_cdr3_a  = 12,
    int weight_cdr3_b       = 3,
    int gap_penalty_cdr3_b  = 12
) {
    // ---- prepare inputs -------------------------------------------------------
    VDistLookup vla, vlb;
    vla.build(v_dist_a);
    vlb.build(v_dist_b);

    PreparedTCRs t = prepare_tcrs(va_genes, cdr3a_seqs, vb_genes, cdr3b_seqs,
                                  vla, vlb, "rcpp_tcrdist_matrix");
    const int n = t.n;

    // ---- allocate result (zero-initialised: diagonal is 0) -----------------
    NumericMatrix result(n, n);

    // ---- upper triangle + mirror -------------------------------------------
    for (int i = 0; i < n - 1; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int ria = t.vi_a[i];
        const int rib = t.vi_b[i];

        for (int j = i + 1; j < n; ++j) {
            const double vd_a = vla.lookup(ria, t.vi_a[j]);
            const double vd_b = vlb.lookup(rib, t.vi_b[j]);

            const double cd_a = cdr3_dist_fast(t.cdr3a[i], t.cdr3a[j],
                                                weight_cdr3_a,
                                                gap_penalty_cdr3_a);
            const double cd_b = cdr3_dist_fast(t.cdr3b[i], t.cdr3b[j],
                                                weight_cdr3_b,
                                                gap_penalty_cdr3_b);

            const double d = vd_a + cd_a + vd_b + cd_b;
            result(i, j) = d;
            result(j, i) = d;
        }
    }

    return result;
}
