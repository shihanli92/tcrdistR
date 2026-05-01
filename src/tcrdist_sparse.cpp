// tcrdist_sparse.cpp — Sparse TCRdist matrix with early termination.
//
// Exports:
//   rcpp_tcrdist_sparse() — COO triplet list of pairs within a distance threshold.
//
// Only pairs (i, j) with i < j and distance <= threshold are stored.
// Three-stage early termination:
//   1. V-dist sum > threshold  -> skip
//   2. V-dist + CDR3-alpha > threshold -> skip
//   3. Full distance > threshold -> skip

#include "tcrdist_core.h"

using namespace Rcpp;

//' Sparse TCRdist matrix as COO triplets (C++ implementation)
//'
//' Computes pairwise TCRdist distances for all (i, j) pairs with i < j and
//' distance <= \code{threshold}, returning a sparse representation in
//' coordinate (COO) format.
//'
//' Three-stage early termination is used to skip pairs that cannot satisfy
//' the threshold:
//' \enumerate{
//'   \item If V-region distance alone exceeds threshold, skip.
//'   \item If V-region + CDR3-alpha distance exceeds threshold, skip.
//'   \item If full distance exceeds threshold, skip.
//' }
//'
//' \code{Rcpp::checkUserInterrupt()} is called every 100 outer-loop rows.
//'
//' @param va_genes   Character vector of length N. Alpha V-gene allele names.
//' @param cdr3a_seqs Character vector of length N. Alpha CDR3 sequences.
//' @param vb_genes   Character vector of length N. Beta V-gene allele names.
//' @param cdr3b_seqs Character vector of length N. Beta CDR3 sequences.
//' @param v_dist_a   Named square \code{NumericMatrix}. Pre-computed pairwise
//'   V-region distances for alpha genes.
//' @param v_dist_b   Named square \code{NumericMatrix}. Pre-computed pairwise
//'   V-region distances for beta genes.
//' @param threshold  Numeric. Maximum distance to include in the result.
//'   Pairs with distance > threshold are omitted.
//' @param weight_cdr3_a Integer. Alpha CDR3 alignment distance multiplier.
//'   Default 3.
//' @param gap_penalty_cdr3_a Integer. Alpha per-residue length-difference
//'   penalty. Default 12.
//' @param weight_cdr3_b Integer. Beta CDR3 alignment distance multiplier.
//'   Default 3.
//' @param gap_penalty_cdr3_b Integer. Beta per-residue length-difference
//'   penalty. Default 12.
//' @return A \code{List} with four elements:
//'   \describe{
//'     \item{\code{i}}{Integer vector of 1-based row indices.}
//'     \item{\code{j}}{Integer vector of 1-based column indices.}
//'     \item{\code{x}}{Numeric vector of distance values.}
//'     \item{\code{n}}{Integer scalar: matrix dimension (N).}
//'   }
//' @examples
//' \dontrun{
//'   triplets <- rcpp_tcrdist_sparse(
//'     va_genes = c("TRAV1-1*01", "TRAV1-2*01"),
//'     cdr3a_seqs = c("CAVSANSGTYF", "CAVSANSGTYF"),
//'     vb_genes = c("TRBV20-1*01", "TRBV20-1*01"),
//'     cdr3b_seqs = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
//'     v_dist_a = v_alpha_mat, v_dist_b = v_beta_mat,
//'     threshold = 50
//'   )
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_tcrdist_sparse(
    const CharacterVector& va_genes,
    const CharacterVector& cdr3a_seqs,
    const CharacterVector& vb_genes,
    const CharacterVector& cdr3b_seqs,
    const NumericMatrix&   v_dist_a,
    const NumericMatrix&   v_dist_b,
    double threshold,
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
                                  vla, vlb, "rcpp_tcrdist_sparse");
    const int n = t.n;

    // ---- COO triplet storage (growing vectors) ------------------------------
    std::vector<int>    row_idx;
    std::vector<int>    col_idx;
    std::vector<double> dist_vals;

    // ---- upper triangle with early termination -----------------------------
    for (int i = 0; i < n - 1; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int ria = t.vi_a[i];
        const int rib = t.vi_b[i];

        for (int j = i + 1; j < n; ++j) {
            // Stage 1: V-region sum only
            const double vd = vla.lookup(ria, t.vi_a[j])
                            + vlb.lookup(rib, t.vi_b[j]);
            if (vd > threshold) continue;

            // Stage 2: add CDR3-alpha
            const double cd_a = cdr3_dist_fast(t.cdr3a[i], t.cdr3a[j],
                                               weight_cdr3_a,
                                               gap_penalty_cdr3_a);
            if (vd + cd_a > threshold) continue;

            // Stage 3: add CDR3-beta for full distance
            const double cd_b = cdr3_dist_fast(t.cdr3b[i], t.cdr3b[j],
                                               weight_cdr3_b,
                                               gap_penalty_cdr3_b);
            const double d = vd + cd_a + cd_b;
            if (d <= threshold) {
                row_idx.push_back(i + 1);   // 1-based
                col_idx.push_back(j + 1);   // 1-based
                dist_vals.push_back(d);
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("i") = Rcpp::wrap(row_idx),
        Rcpp::Named("j") = Rcpp::wrap(col_idx),
        Rcpp::Named("x") = Rcpp::wrap(dist_vals),
        Rcpp::Named("n") = n
    );
}
