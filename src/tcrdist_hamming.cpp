// tcrdist_hamming.cpp — Hamming distance for CDR3 amino acid sequences.
//
// Simple mismatch counting for equal-length sequences.

#include <Rcpp.h>
#include <string>
#include <vector>

// [[Rcpp::export]]
int rcpp_hamming_distance(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) {
        return -1;
    }
    int dist = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            ++dist;
        }
    }
    return dist;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix rcpp_hamming_matrix(const Rcpp::CharacterVector& seqs) {
    int n = seqs.size();
    Rcpp::IntegerMatrix mat(n, n);

    // Precompute std::string copies
    std::vector<std::string> ss(n);
    for (int i = 0; i < n; ++i) {
        ss[i] = Rcpp::as<std::string>(seqs[i]);
    }

    for (int i = 0; i < n; ++i) {
        mat(i, i) = 0;
        for (int j = i + 1; j < n; ++j) {
            int d;
            if (ss[i].size() != ss[j].size()) {
                // Max penalty: sum of both lengths (like a full mismatch)
                d = static_cast<int>(std::max(ss[i].size(), ss[j].size()));
            } else {
                d = 0;
                for (size_t k = 0; k < ss[i].size(); ++k) {
                    if (ss[i][k] != ss[j][k]) {
                        ++d;
                    }
                }
            }
            mat(i, j) = d;
            mat(j, i) = d;
        }
    }

    return mat;
}
