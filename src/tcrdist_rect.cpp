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
//' @export
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
    int weight_cdr3_region       = 3,
    int gap_penalty_cdr3_region  = 12
) {
    const int nq = query_va.size();
    const int nr = ref_va.size();

    // ---- validate input lengths --------------------------------------------
    if (query_cdr3a.size() != nq || query_vb.size() != nq || query_cdr3b.size() != nq) {
        Rcpp::stop(
            "rcpp_tcrdist_rect: query_va, query_cdr3a, query_vb, query_cdr3b "
            "must all have the same length"
        );
    }
    if (ref_cdr3a.size() != nr || ref_vb.size() != nr || ref_cdr3b.size() != nr) {
        Rcpp::stop(
            "rcpp_tcrdist_rect: ref_va, ref_cdr3a, ref_vb, ref_cdr3b "
            "must all have the same length"
        );
    }

    // ---- build V-gene lookup tables ----------------------------------------
    VDistLookup vla, vlb;
    vla.build(v_dist_a);
    vlb.build(v_dist_b);

    // ---- pre-convert CharacterVectors to std::string -----------------------
    std::vector<std::string> qva(nq), qcdr3a(nq), qvb(nq), qcdr3b(nq);
    for (int k = 0; k < nq; ++k) {
        qva[k]    = Rcpp::as<std::string>(query_va[k]);
        qcdr3a[k] = Rcpp::as<std::string>(query_cdr3a[k]);
        qvb[k]    = Rcpp::as<std::string>(query_vb[k]);
        qcdr3b[k] = Rcpp::as<std::string>(query_cdr3b[k]);
    }

    std::vector<std::string> rva(nr), rcdr3a(nr), rvb(nr), rcdr3b(nr);
    for (int k = 0; k < nr; ++k) {
        rva[k]    = Rcpp::as<std::string>(ref_va[k]);
        rcdr3a[k] = Rcpp::as<std::string>(ref_cdr3a[k]);
        rvb[k]    = Rcpp::as<std::string>(ref_vb[k]);
        rcdr3b[k] = Rcpp::as<std::string>(ref_cdr3b[k]);
    }

    // ---- pre-resolve V-gene names to indices -------------------------------
    std::vector<int> qri_a(nq), qri_b(nq);
    for (int k = 0; k < nq; ++k) {
        int ia = vla.resolve(qva[k]);
        if (ia < 0) {
            Rcpp::stop(
                "rcpp_tcrdist_rect: query alpha V-gene '%s' not found in v_dist_a",
                qva[k].c_str()
            );
        }
        qri_a[k] = ia;

        int ib = vlb.resolve(qvb[k]);
        if (ib < 0) {
            Rcpp::stop(
                "rcpp_tcrdist_rect: query beta V-gene '%s' not found in v_dist_b",
                qvb[k].c_str()
            );
        }
        qri_b[k] = ib;
    }

    std::vector<int> rri_a(nr), rri_b(nr);
    for (int k = 0; k < nr; ++k) {
        int ia = vla.resolve(rva[k]);
        if (ia < 0) {
            Rcpp::stop(
                "rcpp_tcrdist_rect: reference alpha V-gene '%s' not found in v_dist_a",
                rva[k].c_str()
            );
        }
        rri_a[k] = ia;

        int ib = vlb.resolve(rvb[k]);
        if (ib < 0) {
            Rcpp::stop(
                "rcpp_tcrdist_rect: reference beta V-gene '%s' not found in v_dist_b",
                rvb[k].c_str()
            );
        }
        rri_b[k] = ib;
    }

    // ---- preprocess all CDR3 sequences -------------------------------------
    std::vector<CDR3Data> qcdr3a_d(nq), qcdr3b_d(nq);
    for (int k = 0; k < nq; ++k) {
        qcdr3a_d[k] = preprocess_cdr3(qcdr3a[k]);
        qcdr3b_d[k] = preprocess_cdr3(qcdr3b[k]);
    }

    std::vector<CDR3Data> rcdr3a_d(nr), rcdr3b_d(nr);
    for (int k = 0; k < nr; ++k) {
        rcdr3a_d[k] = preprocess_cdr3(rcdr3a[k]);
        rcdr3b_d[k] = preprocess_cdr3(rcdr3b[k]);
    }

    // ---- allocate result matrix --------------------------------------------
    NumericMatrix result(nq, nr);

    // ---- full rectangular loop: no symmetry --------------------------------
    for (int i = 0; i < nq; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int qria = qri_a[i];
        const int qrib = qri_b[i];

        for (int j = 0; j < nr; ++j) {
            const double vd_a = vla.lookup(qria, rri_a[j]);
            const double vd_b = vlb.lookup(qrib, rri_b[j]);

            const double cd_a = cdr3_dist_fast(qcdr3a_d[i], rcdr3a_d[j],
                                               weight_cdr3_region,
                                               gap_penalty_cdr3_region);
            const double cd_b = cdr3_dist_fast(qcdr3b_d[i], rcdr3b_d[j],
                                               weight_cdr3_region,
                                               gap_penalty_cdr3_region);

            result(i, j) = vd_a + cd_a + vd_b + cd_b;
        }
    }

    return result;
}
