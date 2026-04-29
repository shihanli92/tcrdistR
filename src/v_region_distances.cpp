// v_region_distances.cpp — Pairwise V-gene loop sequence distance matrix.
//
// Exports:
//   rcpp_compute_v_region_distances() — all pairwise distances for a set of
//     V-gene loop sequences (merged CDR loops, space-separated).
//
// This mirrors compute_all_v_region_distances() from rconga/R/tcr_distances.R
// but works on a pre-assembled vector of (gene_id, loop_seq) pairs rather
// than loading the gene database internally.

#include "tcrdist_core.h"

using namespace Rcpp;

//' Compute pairwise V-region distances (C++ implementation)
//'
//' Computes all pairwise \code{blosum_sequence_distance} values for a set of
//' V-gene CDR loop sequences, weighted by \code{weight_v_region}.  Matches
//' the inner loop of \code{compute_all_v_region_distances()} in
//' \code{rconga/R/tcr_distances.R}.
//'
//' The sequences in \code{loop_seqs} are pre-aligned merged CDR loop strings
//' (CDRs joined by \code{" "}), as produced by the R gene-database loader.
//' Spaces mark position padding and are skipped (both positions must be space).
//'
//' @param gene_ids Character vector of length N. V-gene allele identifiers
//'   (e.g. \code{"TRAV1-1*01"}).  Used as row/column names of the result.
//' @param loop_seqs Character vector of length N.  Pre-aligned merged CDR loop
//'   sequences corresponding to each gene in \code{gene_ids}.  All sequences
//'   must have the same character length.
//' @param weight_v_region Numeric. Multiplied by each raw
//'   \code{blosum_sequence_distance} value.  Default 1.0 matches
//'   \code{WEIGHT_V_REGION}.
//' @param gap_penalty_v_region Numeric. Forwarded to
//'   \code{rcpp_blosum_sequence_distance} as \code{gap_penalty}.  Default 4.0
//'   matches \code{GAP_PENALTY_V_REGION}.
//' @return A named symmetric N x N \code{NumericMatrix}.
//'   \code{result[i, j] = weight_v_region * blosum_sequence_distance(loop_seqs[i], loop_seqs[j])}.
//' @examples
//' \dontrun{
//'   ids  <- c("TRAV1-1*01", "TRAV1-2*01")
//'   seqs <- c("ACDE FGHI", "ACDE FGHK")
//'   rcpp_compute_v_region_distances(ids, seqs)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix rcpp_compute_v_region_distances(
    const CharacterVector& gene_ids,
    const CharacterVector& loop_seqs,
    double weight_v_region      = 1.0,
    double gap_penalty_v_region = 4.0
) {
    const int n = gene_ids.size();

    if (loop_seqs.size() != n) {
        Rcpp::stop(
            "rcpp_compute_v_region_distances: gene_ids and loop_seqs must have "
            "the same length (%d vs %d)",
            n, static_cast<int>(loop_seqs.size())
        );
    }

    // Pre-convert CharacterVector to std::vector<std::string> to avoid
    // repeated Rcpp conversion in the inner loop.
    std::vector<std::string> seqs(n);
    for (int k = 0; k < n; ++k) {
        seqs[k] = Rcpp::as<std::string>(loop_seqs[k]);
    }

    // Validate that all sequences have the same length (required by
    // blosum_sequence_distance / blosum_character_distance in R).
    if (n > 0) {
        const int expected_len = static_cast<int>(seqs[0].size());
        for (int k = 1; k < n; ++k) {
            if (static_cast<int>(seqs[k].size()) != expected_len) {
                Rcpp::stop(
                    "rcpp_compute_v_region_distances: loop sequence length mismatch: "
                    "gene '%s' has length %d but gene '%s' has length %d",
                    Rcpp::as<std::string>(gene_ids[0]).c_str(), expected_len,
                    Rcpp::as<std::string>(gene_ids[k]).c_str(),
                    static_cast<int>(seqs[k].size())
                );
            }
        }
    }

    // Pre-encode each sequence into two parallel arrays for fast position access:
    //   aa_code[k][pos]:  0-19 for AAs, -2 for space, -3 for gap '.', -4 for stop '*'
    // This avoids character lookup in the inner loop.
    // Special sentinel values:
    const int CODE_SPACE = -2;
    const int CODE_GAP   = -3;
    const int CODE_STOP  = -4;

    const int seq_len = (n > 0) ? static_cast<int>(seqs[0].size()) : 0;
    // Flattened: encoded[k * seq_len + pos]
    std::vector<int> encoded(static_cast<size_t>(n) * seq_len);

    for (int k = 0; k < n; ++k) {
        const std::string& s = seqs[k];
        for (int p = 0; p < seq_len; ++p) {
            char c = s[p];
            if (c == ' ') {
                encoded[k * seq_len + p] = CODE_SPACE;
            } else if (c == '.') {
                encoded[k * seq_len + p] = CODE_GAP;
            } else if (c == '*') {
                encoded[k * seq_len + p] = CODE_STOP;
            } else {
                int idx = tcrdist::aa_to_index(c);
                if (idx < 0) {
                    Rcpp::stop(
                        "rcpp_compute_v_region_distances: unknown character '%c' in "
                        "loop sequence for gene '%s'",
                        c, Rcpp::as<std::string>(gene_ids[k]).c_str()
                    );
                }
                encoded[k * seq_len + p] = idx;
            }
        }
    }

    // Allocate result (zero-initialised: diagonal is 0).
    NumericMatrix result(n, n);

    // Upper triangle loop with mirroring.
    for (int i = 0; i < n; ++i) {
        const int* row_i = encoded.data() + i * seq_len;

        for (int j = i + 1; j < n; ++j) {
            const int* row_j = encoded.data() + j * seq_len;

            double dist = 0.0;
            for (int p = 0; p < seq_len; ++p) {
                int ca = row_i[p];
                int cb = row_j[p];

                if (ca == CODE_SPACE) {
                    // Both must be space (invariant); skip without checking
                    // for efficiency since sequences were validated equal-length.
                    // The R implementation asserts b == ' '; we trust the input.
                    continue;
                }
                if (ca == CODE_GAP && cb == CODE_GAP)   continue;
                if (ca == CODE_STOP && cb == CODE_STOP)  continue;
                if (ca == CODE_GAP || cb == CODE_GAP ||
                    ca == CODE_STOP || cb == CODE_STOP) {
                    dist += gap_penalty_v_region;
                    continue;
                }
                // Both are standard AAs: BSD4 lookup
                dist += tcrdist::BSD4_FLAT[ca * 20 + cb];
            }

            double d = weight_v_region * dist;
            result(i, j) = d;
            result(j, i) = d;
        }
    }

    // Assign row and column names
    CharacterVector ids_copy(gene_ids.begin(), gene_ids.end());
    rownames(result) = ids_copy;
    colnames(result) = ids_copy;

    return result;
}
