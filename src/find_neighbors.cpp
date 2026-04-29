// find_neighbors.cpp — KNN from precomputed distance or PCA matrix.
//
// Exports:
//   rcpp_knn_from_distance_matrix() — KNN from N x N distance matrix.
//   rcpp_knn_from_pca_matrix()       — KNN from N x D PCA embeddings.
//
// Adapted from rconga/src/find_neighbors_rcpp.cpp.
// Key changes from rconga:
//   - Output indices are 1-based (R convention).
//   - Named output elements use "knn_indices" / "knn_distances".

#include "tcrdist_core.h"

#include <cmath>

using namespace Rcpp;

//' K-nearest-neighbors from a precomputed distance matrix (C++ implementation)
//'
//' For each row of a square distance matrix \code{D}, finds the K nearest
//' neighbors while masking out entries that share the same alpha-chain or
//' beta-chain group assignment.  Same-group pairs are assigned the sentinel
//' distance 1e3 before selection.
//'
//' Uses \code{std::nth_element} for O(N) partial sorting per row, optionally
//' followed by \code{std::sort} of the K selected neighbors.
//'
//' @param D Numeric matrix (N x N). Precomputed pairwise distance matrix.
//'   Must be square.
//' @param K Integer. Number of nearest neighbors to extract per row.
//' @param agroups Integer vector of length N. Alpha-chain group assignments.
//'   Rows sharing the same agroups value are masked from each other.
//' @param bgroups Integer vector of length N. Beta-chain group assignments.
//'   Rows sharing the same bgroups value are masked from each other.
//' @param sort_nbrs Logical. If \code{TRUE}, sort the K neighbors by ascending
//'   distance.  Default \code{TRUE}.
//' @return A \code{List} with two elements:
//'   \describe{
//'     \item{\code{knn_indices}}{Integer matrix (N x K). 1-based neighbor indices.}
//'     \item{\code{knn_distances}}{Numeric matrix (N x K). Corresponding distances.}
//'   }
//' @examples
//' \dontrun{
//'   D <- matrix(c(0,1,2,1,0,3,2,3,0), nrow=3)
//'   result <- rcpp_knn_from_distance_matrix(D, K=1L, agroups=1:3, bgroups=1:3)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_knn_from_distance_matrix(
    const NumericMatrix& D,
    int K,
    const IntegerVector& agroups,
    const IntegerVector& bgroups,
    bool sort_nbrs = true
) {
    const int N = D.nrow();

    // ---- input validation --------------------------------------------------
    if (D.ncol() != N) {
        Rcpp::stop(
            "rcpp_knn_from_distance_matrix: D must be square (got %d x %d)",
            N, D.ncol()
        );
    }
    if (agroups.size() != N) {
        Rcpp::stop(
            "rcpp_knn_from_distance_matrix: agroups length (%d) != D rows (%d)",
            (int)agroups.size(), N
        );
    }
    if (bgroups.size() != N) {
        Rcpp::stop(
            "rcpp_knn_from_distance_matrix: bgroups length (%d) != D rows (%d)",
            (int)bgroups.size(), N
        );
    }
    if (K <= 0) {
        Rcpp::stop(
            "rcpp_knn_from_distance_matrix: K must be positive (got %d)", K
        );
    }
    if (K >= N) {
        Rcpp::stop(
            "rcpp_knn_from_distance_matrix: K (%d) must be < N (%d)", K, N
        );
    }

    // ---- allocate output matrices ------------------------------------------
    IntegerMatrix nbr_indices(N, K);
    NumericMatrix nbr_distances(N, K);

    static const double MASK_DIST = 1e3;

    // Comparator: ascending distance, tie-break by index
    auto cmp = [](const std::pair<double, int>& a,
                  const std::pair<double, int>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    };

    // Reusable per-row candidate buffer
    std::vector<std::pair<double, int>> candidates(N);

    // ---- per-row KNN extraction --------------------------------------------
    for (int i = 0; i < N; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int ag_i = agroups[i];
        const int bg_i = bgroups[i];

        for (int j = 0; j < N; ++j) {
            double d = D(i, j);
            // Mask same-group (includes self since self shares group with itself)
            if (agroups[j] == ag_i || bgroups[j] == bg_i) {
                d = MASK_DIST;
            }
            candidates[j] = {d, j};
        }

        std::nth_element(candidates.begin(), candidates.begin() + K,
                         candidates.end(), cmp);

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


//' K-nearest-neighbors from a PCA embedding matrix (C++ implementation)
//'
//' Computes K nearest neighbors directly from PCA embeddings (N x D matrix)
//' without materializing the full N x N Euclidean distance matrix.  For each
//' row, computes Euclidean distances to all N points on the fly, applies
//' group masking, and extracts the K nearest via \code{std::nth_element}.
//'
//' Memory: O(N*D + N*K) instead of O(N^2).
//'
//' @param pca_matrix Numeric matrix (N x D). PCA or other embedding coordinates.
//'   Rows are samples, columns are dimensions.
//' @param K Integer. Number of nearest neighbors to extract per point.
//' @param agroups Integer vector of length N. Alpha-chain group assignments.
//' @param bgroups Integer vector of length N. Beta-chain group assignments.
//' @param sort_nbrs Logical. If \code{TRUE}, sort the K neighbors by ascending
//'   Euclidean distance.  Default \code{TRUE}.
//' @return A \code{List} with two elements:
//'   \describe{
//'     \item{\code{knn_indices}}{Integer matrix (N x K). 1-based neighbor indices.}
//'     \item{\code{knn_distances}}{Numeric matrix (N x K). Euclidean distances.}
//'   }
//' @examples
//' \dontrun{
//'   pca <- matrix(rnorm(30), nrow=10, ncol=3)
//'   result <- rcpp_knn_from_pca_matrix(pca, K=3L, agroups=1:10, bgroups=1:10)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rcpp_knn_from_pca_matrix(
    const NumericMatrix& pca_matrix,
    int K,
    const IntegerVector& agroups,
    const IntegerVector& bgroups,
    bool sort_nbrs = true
) {
    const int N = pca_matrix.nrow();
    const int D = pca_matrix.ncol();

    // ---- input validation --------------------------------------------------
    if (agroups.size() != N) {
        Rcpp::stop(
            "rcpp_knn_from_pca_matrix: agroups length (%d) != N (%d)",
            (int)agroups.size(), N
        );
    }
    if (bgroups.size() != N) {
        Rcpp::stop(
            "rcpp_knn_from_pca_matrix: bgroups length (%d) != N (%d)",
            (int)bgroups.size(), N
        );
    }
    if (K <= 0) {
        Rcpp::stop(
            "rcpp_knn_from_pca_matrix: K must be positive (got %d)", K
        );
    }
    if (K >= N) {
        Rcpp::stop(
            "rcpp_knn_from_pca_matrix: K (%d) must be < N (%d)", K, N
        );
    }

    // ---- allocate output matrices ------------------------------------------
    IntegerMatrix nbr_indices(N, K);
    NumericMatrix nbr_distances(N, K);

    static const double MASK_DIST = 1e3;

    // Comparator: ascending distance, tie-break by index
    auto cmp = [](const std::pair<double, int>& a,
                  const std::pair<double, int>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    };

    // Reusable per-row candidate buffer
    std::vector<std::pair<double, int>> candidates(N);

    // Raw pointer to column-major PCA data
    const double* pca_ptr = pca_matrix.begin();

    // ---- per-row KNN computation -------------------------------------------
    for (int i = 0; i < N; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        const int ag_i = agroups[i];
        const int bg_i = bgroups[i];

        for (int j = 0; j < N; ++j) {
            if (agroups[j] == ag_i || bgroups[j] == bg_i) {
                candidates[j] = {MASK_DIST, j};
                continue;
            }
            // Euclidean distance (column-major: element (row, col) = col*N + row)
            double dist_sq = 0.0;
            for (int d = 0; d < D; ++d) {
                double diff = pca_ptr[d * N + i] - pca_ptr[d * N + j];
                dist_sq += diff * diff;
            }
            candidates[j] = {std::sqrt(dist_sq), j};
        }

        std::nth_element(candidates.begin(), candidates.begin() + K,
                         candidates.end(), cmp);

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
