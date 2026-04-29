// blosum.cpp — Rcpp export: expose the constexpr BSD4 table to R.
//
// rcpp_build_bsd4() constructs a named 20x20 NumericMatrix from the
// compile-time tcrdist::BSD4 array so R code can inspect or verify
// the matrix without reimplementing the derivation in R.

#include "tcrdist_core.h"

using namespace Rcpp;

//' Build the BSD4 amino acid distance matrix (C++ constexpr values)
//'
//' Returns the 20x20 BSD4 distance matrix derived from BLOSUM62 using the
//' CoNGA transformation:
//' \itemize{
//'   \item diagonal (a == b): 0
//'   \item off-diagonal, BLOSUM62 < 0: 4
//'   \item off-diagonal, BLOSUM62 >= 0: 4 - BLOSUM62
//' }
//' Row and column names follow alphabetical single-letter amino acid order
//' (A C D E F G H I K L M N P Q R S T V W Y), matching \code{AMINO_ACIDS}
//' in rconga.
//'
//' The values are compiled into the package as \code{constexpr} C++ arrays;
//' this function exposes them to R for verification and for use by any R
//' code that needs the matrix.
//'
//' @return A named 20x20 \code{NumericMatrix} of BSD4 distances.
//' @examples
//' \dontrun{
//'   bsd4 <- rcpp_build_bsd4()
//'   bsd4["A", "S"]  # 3 (BLOSUM62[A,S]=1, BSD4=4-1=3)
//'   bsd4["A", "A"]  # 0 (diagonal)
//'   bsd4["A", "C"]  # 4 (BLOSUM62=-something negative)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix rcpp_build_bsd4() {
    const int n = tcrdist::AA_COUNT;  // 20

    NumericMatrix mat(n, n);

    // Fill from the constexpr BSD4 array (row-major)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mat(i, j) = tcrdist::BSD4[i][j];
        }
    }

    // Assign row and column names from AA_ORDER
    CharacterVector aa_names(n);
    for (int k = 0; k < n; ++k) {
        aa_names[k] = std::string(1, tcrdist::AA_ORDER[k]);
    }

    rownames(mat) = aa_names;
    colnames(mat) = aa_names;

    return mat;
}


//' Look up a single BLOSUM62 score for two amino acids
//'
//' @param aa1 Single-letter amino acid code (length-1 character string).
//' @param aa2 Single-letter amino acid code (length-1 character string).
//' @return Integer BLOSUM62 score for the pair.
//' @keywords internal
// [[Rcpp::export]]
int rcpp_blosum62_lookup(const std::string& aa1, const std::string& aa2) {
    if (aa1.size() != 1u || aa2.size() != 1u) {
        Rcpp::stop("rcpp_blosum62_lookup: both inputs must be single characters");
    }
    int i = tcrdist::aa_to_index(aa1[0]);
    int j = tcrdist::aa_to_index(aa2[0]);
    if (i < 0 || j < 0) {
        Rcpp::stop("rcpp_blosum62_lookup: unknown amino acid '%s' or '%s'",
                    aa1.c_str(), aa2.c_str());
    }
    return tcrdist::BLOSUM62[i][j];
}
