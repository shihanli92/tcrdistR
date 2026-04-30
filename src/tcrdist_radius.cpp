// tcrdist_radius.cpp — Radius-based neighbor search via TCRdist.
//
// Exports:
//   rcpp_tcrdist_radius_neighbors() — variable-length neighbor lists within radius.
//
// Adapted from rconga/src/tcrdist_rcpp.cpp lines 515-628.
// Key changes from rconga:
//   - No bsd4 matrix parameter; uses cdr3_dist_fast() from tcrdist_core.h.
//   - Output indices are 1-based (R convention).
//   - Three-stage early termination: V-dist -> V+CDR3a -> full.
//   - Uses CDR3Data precomputation and VDistLookup from tcrdist_core.h.

#include "tcrdist_core.h"

using namespace Rcpp;

//' Radius-based TCRdist neighbor search with group masking (C++ implementation)
//'
//' For each of the N input TCRs, finds all other TCRs within \code{radius}
//' TCRdist distance, excluding TCRs that share the same alpha or beta group
//' (same-group exclusion) and self-pairs.
//'
//' Three-stage early termination is applied per pair:
//' \enumerate{
//'   \item If V-region distance alone exceeds \code{radius}, skip.
//'   \item If V-region + CDR3-alpha distance exceeds \code{radius}, skip.
//'   \item If full distance exceeds \code{radius}, skip.
//' }
//'
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
//' @param radius     Numeric. Search radius. Only neighbors within this
//'   distance (inclusive) are returned.
//' @param agroups Integer vector of length N. Alpha-chain group assignments.
//' @param bgroups Integer vector of length N. Beta-chain group assignments.
//' @param weight_cdr3_region Integer. CDR3 alignment distance multiplier.
//'   Default 3 matches \code{WEIGHT_CDR3_REGION}.
//' @param gap_penalty_cdr3_region Integer. Per-residue length-difference
//'   penalty. Default 12 matches \code{GAP_PENALTY_CDR3_REGION}.
//' @return A \code{List} of length N. Each element is a \code{List} with:
//'   \describe{
//'     \item{\code{indices}}{Integer vector of 1-based neighbor indices.}
//'     \item{\code{distances}}{Numeric vector of corresponding distances.}
//'   }
//' @examples
//' \dontrun{
//'   result <- rcpp_tcrdist_radius_neighbors(
//'     va_genes = c("TRAV1-1*01", "TRAV1-2*01"),
//'     cdr3a_seqs = c("CAVSANSGTYF", "CAVSANSGTYF"),
//'     vb_genes = c("TRBV20-1*01", "TRBV20-1*01"),
//'     cdr3b_seqs = c("CASSIRSSYEQYF", "CASSIRSYEQYF"),
//'     v_dist_a = v_alpha_mat, v_dist_b = v_beta_mat,
//'     radius = 50, agroups = 1:2, bgroups = 1:2
//'   )
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_tcrdist_radius_neighbors(
    const CharacterVector& va_genes,
    const CharacterVector& cdr3a_seqs,
    const CharacterVector& vb_genes,
    const CharacterVector& cdr3b_seqs,
    const NumericMatrix&   v_dist_a,
    const NumericMatrix&   v_dist_b,
    double radius,
    const IntegerVector&   agroups,
    const IntegerVector&   bgroups,
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
                                  vla, vlb, "rcpp_tcrdist_radius_neighbors");
    const int n = t.n;

    if (agroups.size() != n) {
        Rcpp::stop("rcpp_tcrdist_radius_neighbors: agroups length (%d) != N (%d)",
                   (int)agroups.size(), n);
    }
    if (bgroups.size() != n) {
        Rcpp::stop("rcpp_tcrdist_radius_neighbors: bgroups length (%d) != N (%d)",
                   (int)bgroups.size(), n);
    }

    // ---- per-row radius neighbor search ------------------------------------
    Rcpp::List result(n);

    for (int i = 0; i < n; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int ag_i = agroups[i];
        const int bg_i = bgroups[i];

        std::vector<int>    nbr_idx;
        std::vector<double> nbr_dist;

        for (int j = 0; j < n; ++j) {
            // Skip self and same-group pairs
            if (i == j) continue;
            if (agroups[j] == ag_i || bgroups[j] == bg_i) continue;

            // Stage 1: V-region sum early termination
            const double vd = vla.lookup(t.vi_a[i], t.vi_a[j])
                            + vlb.lookup(t.vi_b[i], t.vi_b[j]);
            if (vd > radius) continue;

            // Stage 2: V + CDR3-alpha early termination
            const double cd_a = cdr3_dist_fast(t.cdr3a[i], t.cdr3a[j],
                                               weight_cdr3_a,
                                               gap_penalty_cdr3_a);
            if (vd + cd_a > radius) continue;

            // Stage 3: full distance
            const double cd_b = cdr3_dist_fast(t.cdr3b[i], t.cdr3b[j],
                                               weight_cdr3_b,
                                               gap_penalty_cdr3_b);
            const double d = vd + cd_a + cd_b;
            if (d <= radius) {
                nbr_idx.push_back(j + 1);   // 1-based
                nbr_dist.push_back(d);
            }
        }

        result[i] = Rcpp::List::create(
            Rcpp::Named("indices")   = Rcpp::wrap(nbr_idx),
            Rcpp::Named("distances") = Rcpp::wrap(nbr_dist)
        );
    }

    return result;
}
