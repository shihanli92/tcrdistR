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
//' @param weight_cdr3_region Integer. CDR3 alignment distance multiplier.
//'   Default 3 matches \code{WEIGHT_CDR3_REGION}.
//' @param gap_penalty_cdr3_region Integer. Per-residue length-difference
//'   penalty. Default 12 matches \code{GAP_PENALTY_CDR3_REGION}.
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
//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_tcrdist_sparse(
    const CharacterVector& va_genes,
    const CharacterVector& cdr3a_seqs,
    const CharacterVector& vb_genes,
    const CharacterVector& cdr3b_seqs,
    const NumericMatrix&   v_dist_a,
    const NumericMatrix&   v_dist_b,
    double threshold,
    int weight_cdr3_region       = 3,
    int gap_penalty_cdr3_region  = 12
) {
    const int n = va_genes.size();

    // ---- validate input lengths --------------------------------------------
    if (cdr3a_seqs.size() != n || vb_genes.size() != n || cdr3b_seqs.size() != n) {
        Rcpp::stop(
            "rcpp_tcrdist_sparse: va_genes, cdr3a_seqs, vb_genes, cdr3b_seqs "
            "must all have the same length"
        );
    }

    // ---- build V-gene lookup tables ----------------------------------------
    VDistLookup vla, vlb;
    vla.build(v_dist_a);
    vlb.build(v_dist_b);

    // ---- pre-convert CharacterVector to std::string ------------------------
    std::vector<std::string> va(n), cdr3a(n), vb(n), cdr3b(n);
    for (int k = 0; k < n; ++k) {
        va[k]    = Rcpp::as<std::string>(va_genes[k]);
        cdr3a[k] = Rcpp::as<std::string>(cdr3a_seqs[k]);
        vb[k]    = Rcpp::as<std::string>(vb_genes[k]);
        cdr3b[k] = Rcpp::as<std::string>(cdr3b_seqs[k]);
    }

    // ---- pre-resolve V-gene names to indices -------------------------------
    std::vector<int> ri_a(n), ri_b(n);
    for (int k = 0; k < n; ++k) {
        int ia = vla.resolve(va[k]);
        if (ia < 0) {
            Rcpp::stop(
                "rcpp_tcrdist_sparse: alpha V-gene '%s' not found in v_dist_a",
                va[k].c_str()
            );
        }
        ri_a[k] = ia;

        int ib = vlb.resolve(vb[k]);
        if (ib < 0) {
            Rcpp::stop(
                "rcpp_tcrdist_sparse: beta V-gene '%s' not found in v_dist_b",
                vb[k].c_str()
            );
        }
        ri_b[k] = ib;
    }

    // ---- preprocess all CDR3 sequences -------------------------------------
    std::vector<CDR3Data> cdr3a_d(n), cdr3b_d(n);
    for (int k = 0; k < n; ++k) {
        cdr3a_d[k] = preprocess_cdr3(cdr3a[k]);
        cdr3b_d[k] = preprocess_cdr3(cdr3b[k]);
    }

    // ---- COO triplet storage (growing vectors) ------------------------------
    std::vector<int>    row_idx;
    std::vector<int>    col_idx;
    std::vector<double> dist_vals;

    // ---- upper triangle with early termination -----------------------------
    for (int i = 0; i < n - 1; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int ria = ri_a[i];
        const int rib = ri_b[i];

        for (int j = i + 1; j < n; ++j) {
            // Stage 1: V-region sum only
            const double vd = vla.lookup(ria, ri_a[j]) + vlb.lookup(rib, ri_b[j]);
            if (vd > threshold) continue;

            // Stage 2: add CDR3-alpha
            const double cd_a = cdr3_dist_fast(cdr3a_d[i], cdr3a_d[j],
                                               weight_cdr3_region,
                                               gap_penalty_cdr3_region);
            if (vd + cd_a > threshold) continue;

            // Stage 3: add CDR3-beta for full distance
            const double cd_b = cdr3_dist_fast(cdr3b_d[i], cdr3b_d[j],
                                               weight_cdr3_region,
                                               gap_penalty_cdr3_region);
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
