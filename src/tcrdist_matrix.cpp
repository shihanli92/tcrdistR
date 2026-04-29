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
//' @param weight_cdr3_region Integer. CDR3 alignment distance multiplier.
//'   Default 3 matches \code{WEIGHT_CDR3_REGION}.
//' @param gap_penalty_cdr3_region Integer. Per-residue length-difference
//'   penalty.  Default 12 matches \code{GAP_PENALTY_CDR3_REGION}.
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
//' @export
// [[Rcpp::export]]
NumericMatrix rcpp_tcrdist_matrix(
    const CharacterVector& va_genes,
    const CharacterVector& cdr3a_seqs,
    const CharacterVector& vb_genes,
    const CharacterVector& cdr3b_seqs,
    const NumericMatrix&   v_dist_a,
    const NumericMatrix&   v_dist_b,
    int weight_cdr3_region       = 3,
    int gap_penalty_cdr3_region  = 12
) {
    const int n = va_genes.size();

    // ---- validate input lengths --------------------------------------------
    if (cdr3a_seqs.size() != n || vb_genes.size() != n || cdr3b_seqs.size() != n) {
        Rcpp::stop(
            "rcpp_tcrdist_matrix: va_genes, cdr3a_seqs, vb_genes, cdr3b_seqs "
            "must all have the same length"
        );
    }

    // ---- build V-gene lookup tables ----------------------------------------
    VDistLookup vla, vlb;
    vla.build(v_dist_a);
    vlb.build(v_dist_b);

    // ---- pre-convert CharacterVector to std::string for inner loop ---------
    // Avoids repeated Rcpp string conversion overhead.
    std::vector<std::string> va(n), cdr3a(n), vb(n), cdr3b(n);
    for (int k = 0; k < n; ++k) {
        va[k]    = Rcpp::as<std::string>(va_genes[k]);
        cdr3a[k] = Rcpp::as<std::string>(cdr3a_seqs[k]);
        vb[k]    = Rcpp::as<std::string>(vb_genes[k]);
        cdr3b[k] = Rcpp::as<std::string>(cdr3b_seqs[k]);
    }

    // ---- pre-resolve all V-gene names to integer indices -------------------
    // Fails fast with an informative error rather than silently in the inner loop.
    std::vector<int> ri_a(n), ri_b(n);
    for (int k = 0; k < n; ++k) {
        int ia = vla.resolve(va[k]);
        if (ia < 0) {
            Rcpp::stop(
                "rcpp_tcrdist_matrix: alpha V-gene '%s' not found in v_dist_a",
                va[k].c_str()
            );
        }
        ri_a[k] = ia;

        int ib = vlb.resolve(vb[k]);
        if (ib < 0) {
            Rcpp::stop(
                "rcpp_tcrdist_matrix: beta V-gene '%s' not found in v_dist_b",
                vb[k].c_str()
            );
        }
        ri_b[k] = ib;
    }

    // ---- preprocess all CDR3 sequences -------------------------------------
    // Converts AA characters to uint8_t indices and computes gappos/remainder.
    std::vector<CDR3Data> cdr3a_d(n), cdr3b_d(n);
    for (int k = 0; k < n; ++k) {
        cdr3a_d[k] = preprocess_cdr3(cdr3a[k]);
        cdr3b_d[k] = preprocess_cdr3(cdr3b[k]);
    }

    // ---- allocate result (zero-initialised: diagonal is 0) -----------------
    NumericMatrix result(n, n);

    // ---- upper triangle + mirror -------------------------------------------
    for (int i = 0; i < n - 1; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int ria = ri_a[i];
        const int rib = ri_b[i];

        for (int j = i + 1; j < n; ++j) {
            // V-region distances: O(1) flat array lookup
            const double vd_a = vla.lookup(ria, ri_a[j]);
            const double vd_b = vlb.lookup(rib, ri_b[j]);

            // CDR3 distances: uses constexpr BSD4_FLAT
            const double cd_a = cdr3_dist_fast(cdr3a_d[i], cdr3a_d[j],
                                                weight_cdr3_region,
                                                gap_penalty_cdr3_region);
            const double cd_b = cdr3_dist_fast(cdr3b_d[i], cdr3b_d[j],
                                                weight_cdr3_region,
                                                gap_penalty_cdr3_region);

            const double d = vd_a + cd_a + vd_b + cd_b;
            result(i, j) = d;
            result(j, i) = d;
        }
    }

    return result;
}
