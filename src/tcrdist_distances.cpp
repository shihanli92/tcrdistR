// tcrdist_distances.cpp — Single-pair distance functions.
//
// Exports:
//   rcpp_weighted_cdr3_distance()    — CDR3 component distance
//   rcpp_blosum_sequence_distance()  — aligned V-region sequence distance
//
// Both functions produce bit-for-bit identical output to their R equivalents
// in rconga/R/tcr_distances.R.

#include "tcrdist_core.h"

using namespace Rcpp;

// ---------------------------------------------------------------------------
// rcpp_weighted_cdr3_distance
// ---------------------------------------------------------------------------

//' Weighted CDR3 distance (C++ implementation)
//'
//' Direct C++ port of the CDR3 distance algorithm from
//' \code{rconga/src/tcrdist_rcpp.cpp}.  Uses hard-coded \code{ntrim = 3},
//' \code{ctrim = 2} (\code{TRIM_CDR3S = TRUE}) and the fixed gap position
//' formula \code{min(6, 3 + (lenshort - 5) %/% 2)} (\code{ALIGN_CDR3S =
//' FALSE}).
//'
//' BSD4 values are compiled in as \code{constexpr} (no matrix argument
//' needed).  The result is numerically identical to
//' \code{rconga::rcpp_weighted_cdr3_distance()} for valid inputs.
//'
//' @param seq1 Character string. First CDR3 amino acid sequence.
//' @param seq2 Character string. Second CDR3 amino acid sequence.
//' @param weight_cdr3_region Integer. Multiplied by the raw alignment distance.
//'   Default 3 matches \code{WEIGHT_CDR3_REGION}.
//' @param gap_penalty_cdr3_region Integer. Per-residue length-difference penalty.
//'   Default 12 matches \code{GAP_PENALTY_CDR3_REGION}.
//' @return Numeric scalar: the weighted CDR3 distance.
//' @examples
//' \dontrun{
//'   rcpp_weighted_cdr3_distance("CASSIRSSYEQYF", "CASSIRSYEQYF")
//' }
//' @export
// [[Rcpp::export]]
double rcpp_weighted_cdr3_distance(
    const std::string& seq1,
    const std::string& seq2,
    int weight_cdr3_region       = 3,
    int gap_penalty_cdr3_region  = 12
) {
    // Preprocess converts AA chars to indices and computes gappos/remainder.
    // Stops with an informative error for invalid amino acid characters.
    CDR3Data d1 = preprocess_cdr3(seq1);
    CDR3Data d2 = preprocess_cdr3(seq2);

    return cdr3_dist_fast(d1, d2, weight_cdr3_region, gap_penalty_cdr3_region);
}


// ---------------------------------------------------------------------------
// rcpp_blosum_sequence_distance
// ---------------------------------------------------------------------------

//' Aligned sequence distance using BSD4 substitution costs (C++ implementation)
//'
//' Computes the total BSD4 distance between two pre-aligned amino acid strings
//' of equal length.  Matches \code{blosum_sequence_distance()} in
//' \code{rconga/R/tcr_distances.R} exactly.
//'
//' Position handling:
//' \itemize{
//'   \item Both characters are \code{" "} (space): position skipped (padding).
//'     The function stops with an error if only one character is a space.
//'   \item Both characters are \code{"."} (gap): distance 0.
//'   \item Both characters are \code{"*"} (stop codon): distance 0.
//'   \item One character is \code{"."} or \code{"*"}: \code{gap_penalty}.
//'   \item Both are valid amino acids: BSD4 lookup from the constexpr table.
//' }
//'
//' @param seq1 Character string. Aligned sequence (may contain \code{" "},
//'   \code{"."}, or \code{"*"}).
//' @param seq2 Character string. Aligned sequence of the same length as
//'   \code{seq1}.
//' @param gap_penalty Numeric. Applied when exactly one position is a gap or
//'   stop-codon character. Default 4.0 matches \code{GAP_PENALTY_V_REGION}.
//' @return Numeric scalar: sum of per-position BSD4 distances.
//' @examples
//' \dontrun{
//'   rcpp_blosum_sequence_distance("ACDE", "ACDF")
//' }
//' @export
// [[Rcpp::export]]
double rcpp_blosum_sequence_distance(
    const std::string& seq1,
    const std::string& seq2,
    double gap_penalty = 4.0
) {
    if (seq1.size() != seq2.size()) {
        Rcpp::stop(
            "rcpp_blosum_sequence_distance: sequences must have equal length (%d vs %d)",
            static_cast<int>(seq1.size()),
            static_cast<int>(seq2.size())
        );
    }

    const int len = static_cast<int>(seq1.size());
    const char GAP_CHAR = '.';   // GAP_CHARACTER in rconga R
    const char STOP_CHAR = '*';

    double dist = 0.0;

    for (int i = 0; i < len; ++i) {
        char a = seq1[i];
        char b = seq2[i];

        // Padding spaces: both must be space (invariant from Python source)
        if (a == ' ') {
            if (b != ' ') {
                Rcpp::stop(
                    "rcpp_blosum_sequence_distance: space mismatch at position %d "
                    "(seq1=' ', seq2='%c')", i + 1, b
                );
            }
            continue;  // skip position
        }

        // Both gaps -> 0
        if (a == GAP_CHAR && b == GAP_CHAR) continue;
        // Both stop codons -> 0
        if (a == STOP_CHAR && b == STOP_CHAR) continue;

        // One gap or stop codon -> gap penalty
        if (a == GAP_CHAR || b == GAP_CHAR ||
            a == STOP_CHAR || b == STOP_CHAR) {
            dist += gap_penalty;
            continue;
        }

        // Standard amino acid pair: BSD4 lookup via constexpr aa_to_index
        int ia = tcrdist::aa_to_index(a);
        int ib = tcrdist::aa_to_index(b);
        if (ia < 0) {
            Rcpp::stop(
                "rcpp_blosum_sequence_distance: unknown amino acid '%c' at position %d in seq1",
                a, i + 1
            );
        }
        if (ib < 0) {
            Rcpp::stop(
                "rcpp_blosum_sequence_distance: unknown amino acid '%c' at position %d in seq2",
                b, i + 1
            );
        }
        dist += tcrdist::BSD4_FLAT[ia * 20 + ib];
    }

    return dist;
}
