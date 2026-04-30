// tcrdist_rect.cpp — Rectangular query-vs-reference TCRdist matrix.
//
// Exports:
//   rcpp_tcrdist_rect() — nq x nr distance matrix between two TCR sets.
//
// Computes pairwise TCRdist distances between a query set (nq TCRs) and a
// reference set (nr TCRs).  Unlike rcpp_tcrdist_matrix(), no symmetry is
// assumed: all nq * nr pairs are evaluated.

#include "tcrdist_core.h"

using namespace Rcpp;

//' Rectangular TCRdist query-vs-reference distance matrix (C++ implementation)
//'
//' Computes an nq x nr matrix of paired-chain TCRdist distances between a
//' query set of \code{nq} TCRs and a reference set of \code{nr} TCRs.
//' For each query-reference pair \code{(i, j)}:
//' \deqn{
//'   d(i,j) = v\_dist\_a[va\_query_i, va\_ref_j]
//'           + cdr3\_dist(cdr3a\_query_i, cdr3a\_ref_j)
//'           + v\_dist\_b[vb\_query_i, vb\_ref_j]
//'           + cdr3\_dist(cdr3b\_query_i, cdr3b\_ref_j)
//' }
//'
//' No symmetry is exploited; all nq * nr distances are computed.
//' \code{Rcpp::checkUserInterrupt()} is called every 100 outer-loop rows.
//'
//' @param query_va   Character vector of length nq. Query alpha V-gene alleles.
//' @param query_cdr3a Character vector of length nq. Query alpha CDR3 sequences.
//' @param query_vb   Character vector of length nq. Query beta V-gene alleles.
//' @param query_cdr3b Character vector of length nq. Query beta CDR3 sequences.
//' @param ref_va     Character vector of length nr. Reference alpha V-gene alleles.
//' @param ref_cdr3a  Character vector of length nr. Reference alpha CDR3 sequences.
//' @param ref_vb     Character vector of length nr. Reference beta V-gene alleles.
//' @param ref_cdr3b  Character vector of length nr. Reference beta CDR3 sequences.
//' @param v_dist_a   Named square \code{NumericMatrix}. Pre-computed pairwise
//'   V-region distances for alpha genes.
//' @param v_dist_b   Named square \code{NumericMatrix}. Pre-computed pairwise
//'   V-region distances for beta genes.
//' @param weight_cdr3_region Integer. CDR3 alignment distance multiplier.
//'   Default 3 matches \code{WEIGHT_CDR3_REGION}.
//' @param gap_penalty_cdr3_region Integer. Per-residue length-difference
//'   penalty. Default 12 matches \code{GAP_PENALTY_CDR3_REGION}.
//' @return An nq x nr \code{NumericMatrix} of TCRdist distances.
//' @examples
//' \dontrun{
//'   mat <- rcpp_tcrdist_rect(
//'     query_va = c("TRAV1-1*01"),
//'     query_cdr3a = c("CAVSANSGTYF"),
//'     query_vb = c("TRBV20-1*01"),
//'     query_cdr3b = c("CASSIRSSYEQYF"),
//'     ref_va = c("TRAV1-1*01", "TRAV1-2*01"),
//'     ref_cdr3a = c("CAVSANSGTYF", "CAVSANSGTYF"),
//'     ref_vb = c("TRBV20-1*01", "TRBV20-1*01"),
//'     ref_cdr3b = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
//'     v_dist_a = v_alpha_mat,
//'     v_dist_b = v_beta_mat
//'   )
//' }
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix rcpp_tcrdist_rect(
    const CharacterVector& query_va,
    const CharacterVector& query_cdr3a,
    const CharacterVector& query_vb,
    const CharacterVector& query_cdr3b,
    const CharacterVector& ref_va,
    const CharacterVector& ref_cdr3a,
    const CharacterVector& ref_vb,
    const CharacterVector& ref_cdr3b,
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

    PreparedTCRs q = prepare_tcrs(query_va, query_cdr3a, query_vb, query_cdr3b,
                                  vla, vlb, "rcpp_tcrdist_rect (query)");
    PreparedTCRs r = prepare_tcrs(ref_va, ref_cdr3a, ref_vb, ref_cdr3b,
                                  vla, vlb, "rcpp_tcrdist_rect (reference)");

    // ---- allocate result matrix --------------------------------------------
    NumericMatrix result(q.n, r.n);

    // ---- full rectangular loop: no symmetry --------------------------------
    for (int i = 0; i < q.n; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int qria = q.vi_a[i];
        const int qrib = q.vi_b[i];

        for (int j = 0; j < r.n; ++j) {
            const double vd_a = vla.lookup(qria, r.vi_a[j]);
            const double vd_b = vlb.lookup(qrib, r.vi_b[j]);

            const double cd_a = cdr3_dist_fast(q.cdr3a[i], r.cdr3a[j],
                                               weight_cdr3_a,
                                               gap_penalty_cdr3_a);
            const double cd_b = cdr3_dist_fast(q.cdr3b[i], r.cdr3b[j],
                                               weight_cdr3_b,
                                               gap_penalty_cdr3_b);

            result(i, j) = vd_a + cd_a + vd_b + cd_b;
        }
    }

    return result;
}
