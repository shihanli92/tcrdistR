// tcrdist_knn.cpp — K-nearest-neighbors via TCRdist with group masking.
//
// Exports:
//   rcpp_tcrdist_knn() — N x K nearest-neighbor indices and distances.
//
// Adapted from rconga/src/tcrdist_rcpp.cpp lines 377-507.
// Key changes from rconga:
//   - No bsd4 matrix parameter; uses cdr3_dist_fast() from tcrdist_core.h
//     which references the constexpr BSD4_FLAT table.
//   - Output indices are 1-based (R convention).
//   - Uses CDR3Data precomputation and VDistLookup from tcrdist_core.h.

#include "tcrdist_core.h"

using namespace Rcpp;

//' K-nearest-neighbors by TCRdist with group masking (C++ implementation)
//'
//' For each of the N input TCRs, finds the K nearest neighbors by TCRdist
//' while masking out TCRs that share the same alpha-group or beta-group
//' assignment (same-group exclusion).  Same-group pairs and self-pairs are
//' assigned the sentinel distance 1e3 before selection.
//'
//' Uses \code{std::nth_element} for O(N) partial sorting per row, optionally
//' followed by O(K log K) full sorting of the K selected neighbors.
//' \code{Rcpp::checkUserInterrupt()} is called every 100 rows.
//'
//' @param va_genes   Character vector of length N. Alpha V-gene allele names.
//' @param cdr3a_seqs Character vector of length N. Alpha CDR3 sequences.
//' @param vb_genes   Character vector of length N. Beta V-gene allele names.
//' @param cdr3b_seqs Character vector of length N. Beta CDR3 sequences.
//' @param v_dist_a   Named square \code{NumericMatrix}. Pre-computed pairwise
//'   V-region distances for alpha genes.
//' @param v_dist_b   Named square \code{NumericMatrix}. Pre-computed pairwise
//'   V-region distances for beta genes.
//' @param K Integer. Number of nearest neighbors to return per TCR.
//' @param agroups Integer vector of length N. Alpha-chain group assignments.
//'   TCRs sharing the same agroups value are masked from each other.
//' @param bgroups Integer vector of length N. Beta-chain group assignments.
//'   TCRs sharing the same bgroups value are masked from each other.
//' @param sort_nbrs Logical. If \code{TRUE}, sort the K neighbors by ascending
//'   distance. Default \code{TRUE}.
//' @param weight_cdr3_region Integer. CDR3 alignment distance multiplier.
//'   Default 3 matches \code{WEIGHT_CDR3_REGION}.
//' @param gap_penalty_cdr3_region Integer. Per-residue length-difference
//'   penalty. Default 12 matches \code{GAP_PENALTY_CDR3_REGION}.
//' @return A \code{List} with two elements:
//'   \describe{
//'     \item{\code{knn_indices}}{Integer matrix (N x K). 1-based neighbor indices.}
//'     \item{\code{knn_distances}}{Numeric matrix (N x K). Corresponding distances.}
//'   }
//' @examples
//' \dontrun{
//'   result <- rcpp_tcrdist_knn(
//'     va_genes = c("TRAV1-1*01", "TRAV1-2*01"),
//'     cdr3a_seqs = c("CAVSANSGTYF", "CAVSANSGTYF"),
//'     vb_genes = c("TRBV20-1*01", "TRBV20-1*01"),
//'     cdr3b_seqs = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
//'     v_dist_a = v_alpha_mat, v_dist_b = v_beta_mat,
//'     K = 1L,
//'     agroups = 1:2, bgroups = 1:2
//'   )
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_tcrdist_knn(
    const CharacterVector& va_genes,
    const CharacterVector& cdr3a_seqs,
    const CharacterVector& vb_genes,
    const CharacterVector& cdr3b_seqs,
    const NumericMatrix&   v_dist_a,
    const NumericMatrix&   v_dist_b,
    int K,
    const IntegerVector&   agroups,
    const IntegerVector&   bgroups,
    bool sort_nbrs               = true,
    int weight_cdr3_region       = 3,
    int gap_penalty_cdr3_region  = 12
) {
    const int n = va_genes.size();

    // ---- validate inputs ---------------------------------------------------
    if (cdr3a_seqs.size() != n || vb_genes.size() != n || cdr3b_seqs.size() != n) {
        Rcpp::stop(
            "rcpp_tcrdist_knn: va_genes, cdr3a_seqs, vb_genes, cdr3b_seqs "
            "must all have the same length"
        );
    }
    if (agroups.size() != n) {
        Rcpp::stop(
            "rcpp_tcrdist_knn: agroups length (%d) != N (%d)",
            (int)agroups.size(), n
        );
    }
    if (bgroups.size() != n) {
        Rcpp::stop(
            "rcpp_tcrdist_knn: bgroups length (%d) != N (%d)",
            (int)bgroups.size(), n
        );
    }
    if (K <= 0 || K >= n) {
        Rcpp::stop(
            "rcpp_tcrdist_knn: K must be in [1, N-1]; got K=%d, N=%d", K, n
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
                "rcpp_tcrdist_knn: alpha V-gene '%s' not found in v_dist_a",
                va[k].c_str()
            );
        }
        ri_a[k] = ia;

        int ib = vlb.resolve(vb[k]);
        if (ib < 0) {
            Rcpp::stop(
                "rcpp_tcrdist_knn: beta V-gene '%s' not found in v_dist_b",
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

    // ---- allocate output matrices ------------------------------------------
    IntegerMatrix nbr_indices(n, K);
    NumericMatrix nbr_distances(n, K);

    static const double MASK_DIST = 1e3;

    // Comparator: ascending distance, tie-break by index
    auto cmp = [](const std::pair<double, int>& a,
                  const std::pair<double, int>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    };

    // Reusable per-row candidate buffer
    std::vector<std::pair<double, int>> candidates(n);

    // ---- per-row KNN computation -------------------------------------------
    for (int i = 0; i < n; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int ag_i = agroups[i];
        const int bg_i = bgroups[i];

        for (int j = 0; j < n; ++j) {
            if (i == j || agroups[j] == ag_i || bgroups[j] == bg_i) {
                candidates[j] = {MASK_DIST, j};
            } else {
                double d = vla.lookup(ri_a[i], ri_a[j])
                         + cdr3_dist_fast(cdr3a_d[i], cdr3a_d[j],
                                          weight_cdr3_region, gap_penalty_cdr3_region)
                         + vlb.lookup(ri_b[i], ri_b[j])
                         + cdr3_dist_fast(cdr3b_d[i], cdr3b_d[j],
                                          weight_cdr3_region, gap_penalty_cdr3_region);
                candidates[j] = {d, j};
            }
        }

        // Partial sort: K smallest in [0, K)
        std::nth_element(candidates.begin(), candidates.begin() + K,
                         candidates.end(), cmp);

        // Optionally sort the K selected neighbors by distance
        if (sort_nbrs) {
            std::sort(candidates.begin(), candidates.begin() + K, cmp);
        }

        // Store with 1-based indexing
        for (int k = 0; k < K; ++k) {
            nbr_indices(i, k)   = candidates[k].second + 1;  // 1-based
            nbr_distances(i, k) = candidates[k].first;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("knn_indices")   = nbr_indices,
        Rcpp::Named("knn_distances") = nbr_distances
    );
}
