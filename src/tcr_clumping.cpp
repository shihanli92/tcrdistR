// tcrdistR: Poisson clumping test for TCR neighbor enrichment
// Ported from rconga/src/tcr_clumping_rcpp.cpp

#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <string>

// [[Rcpp::export]]
Rcpp::List rcpp_poisson_test_loop(
    const Rcpp::List& all_nbr_indices,
    const Rcpp::List& all_nbr_distances,
    const Rcpp::NumericMatrix& bg_freqs,
    const Rcpp::IntegerVector& agroups,
    const Rcpp::IntegerVector& bgroups,
    const Rcpp::IntegerVector& radii,
    double n_bg_pairs,
    double pvalue_threshold,
    int num_clones,
    Rcpp::Nullable<Rcpp::IntegerVector> clusters_gex_nullable,
    bool use_conservative_pvalues)
{
    int num_radii = radii.size();
    bool has_clusters = clusters_gex_nullable.isNotNull();
    Rcpp::IntegerVector clusters_gex;
    if (has_clusters) {
        clusters_gex = Rcpp::as<Rcpp::IntegerVector>(clusters_gex_nullable);
        if (clusters_gex.size() != num_clones) {
            Rcpp::stop("rcpp_poisson_test_loop: clusters_gex length (%d) != num_clones (%d)",
                        (int)clusters_gex.size(), num_clones);
        }
    }

    // Pre-compute max_nbrs per clone (hoist out of radius loop)
    std::vector<int> max_nbrs(num_clones);
    for (int ii = 0; ii < num_clones; ii++) {
        int count = 0;
        for (int jj = 0; jj < num_clones; jj++) {
            if (jj != ii && agroups[jj] != agroups[ii] &&
                bgroups[jj] != bgroups[ii]) {
                count++;
            }
        }
        max_nbrs[ii] = count;
    }

    // Pre-compute cluster sizes and per-clone intra-cluster max_nbrs
    // max_nbrs_intra[ii] = count of clones in same cluster with different groups
    std::vector<int> max_nbrs_intra;
    std::vector<int> cluster_sizes;
    if (has_clusters) {
        int max_cluster = 0;
        for (int ii = 0; ii < num_clones; ii++) {
            if (clusters_gex[ii] > max_cluster)
                max_cluster = clusters_gex[ii];
        }
        cluster_sizes.resize(max_cluster + 1, 0);
        for (int ii = 0; ii < num_clones; ii++) {
            cluster_sizes[clusters_gex[ii]]++;
        }

        max_nbrs_intra.resize(num_clones);
        for (int ii = 0; ii < num_clones; ii++) {
            int cl = clusters_gex[ii];
            int count = 0;
            for (int jj = 0; jj < num_clones; jj++) {
                if (jj != ii && clusters_gex[jj] == cl &&
                    agroups[jj] != agroups[ii] &&
                    bgroups[jj] != bgroups[ii]) {
                    count++;
                }
            }
            max_nbrs_intra[ii] = count;
        }
    }

    // Result accumulators
    std::vector<int>    res_clone_index;
    std::vector<int>    res_nbr_radius;
    std::vector<double> res_pvalue_adj;
    std::vector<int>    res_num_nbrs;
    std::vector<double> res_expected;
    std::vector<double> res_raw_count;
    std::vector<std::string> res_clump_type;

    Rcpp::LogicalVector is_clumped(num_clones, false);
    Rcpp::NumericMatrix all_raw_pvalues(num_clones, num_radii);
    // Initialize to 1.0
    for (int i = 0; i < num_clones; i++)
        for (int j = 0; j < num_radii; j++)
            all_raw_pvalues(i, j) = 1.0;

    if (all_nbr_indices.size() != num_clones) {
        Rcpp::stop("rcpp_poisson_test_loop: all_nbr_indices length (%d) != num_clones (%d)",
                    (int)all_nbr_indices.size(), num_clones);
    }
    if (all_nbr_distances.size() != num_clones) {
        Rcpp::stop("rcpp_poisson_test_loop: all_nbr_distances length (%d) != num_clones (%d)",
                    (int)all_nbr_distances.size(), num_clones);
    }

    for (int ii = 0; ii < num_clones; ii++) {
        Rcpp::IntegerVector nbr_idx = all_nbr_indices[ii];
        Rcpp::NumericVector nbr_dist = all_nbr_distances[ii];
        int n_nbrs_total = nbr_idx.size();

        for (int irad = 0; irad < num_radii; irad++) {
            int radius = radii[irad];

            // Count neighbors within this radius
            int num_nbrs = 0;
            for (int k = 0; k < n_nbrs_total; k++) {
                if (nbr_dist[k] <= radius) num_nbrs++;
            }
            if (num_nbrs < 1) continue;

            // bg_freqs columns are 0-indexed in C++: column (radius) = P(dist <= radius)
            // But bg_freqs is passed as R matrix (1-indexed columns), so column index
            // for radius d is d (0-based in C++ since NumericMatrix is 0-based)
            double freq = bg_freqs(ii, radius);
            double mu = max_nbrs[ii] * freq;

            // ppois(num_nbrs - 1, mu, lower.tail = FALSE) = P(X >= num_nbrs)
            double pval = R::ppois(num_nbrs - 1, mu, 0, 0);  // lower=FALSE, log=FALSE
            all_raw_pvalues(ii, irad) = pval;

            double pval_adj = pval * num_radii * num_clones;
            if (pval_adj <= pvalue_threshold) {
                is_clumped[ii] = true;
                double raw_count = freq * n_bg_pairs;

                res_clone_index.push_back(ii);  // 0-based
                res_nbr_radius.push_back(radius);
                res_pvalue_adj.push_back(pval_adj);
                res_num_nbrs.push_back(num_nbrs);
                res_expected.push_back(mu);
                res_raw_count.push_back(raw_count);
                res_clump_type.push_back("global");
            }

            // Intra-cluster clumping
            if (has_clusters) {
                int ii_cluster = clusters_gex[ii];
                int num_nbrs_intra = 0;
                for (int k = 0; k < n_nbrs_total; k++) {
                    int nbr = nbr_idx[k];
                    if (nbr < 0 || nbr >= num_clones) {
                        Rcpp::stop("rcpp_poisson_test_loop: invalid neighbor index %d "
                                   "(num_clones=%d) for clone %d", nbr, num_clones, ii);
                    }
                    if (nbr_dist[k] <= radius &&
                        clusters_gex[nbr] == ii_cluster) {
                        num_nbrs_intra++;
                    }
                }
                if (num_nbrs_intra < 1) continue;

                double mu_intra = max_nbrs_intra[ii] * freq;
                double pval_intra;
                if (use_conservative_pvalues) {
                    pval_intra = R::ppois(num_nbrs_intra - 1, mu_intra, 0, 0) *
                                 num_radii * num_clones;
                } else {
                    pval_intra = R::ppois(num_nbrs_intra - 1, mu_intra, 0, 0) *
                                 num_radii * cluster_sizes[ii_cluster];
                }

                if (pval_intra <= pvalue_threshold) {
                    is_clumped[ii] = true;
                    double raw_count_intra = freq * n_bg_pairs;

                    res_clone_index.push_back(ii);
                    res_nbr_radius.push_back(radius);
                    res_pvalue_adj.push_back(pval_intra);
                    res_num_nbrs.push_back(num_nbrs_intra);
                    res_expected.push_back(mu_intra);
                    res_raw_count.push_back(raw_count_intra);
                    res_clump_type.push_back("intra_gex_cluster");
                }
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("clone_index") = Rcpp::wrap(res_clone_index),
        Rcpp::Named("nbr_radius") = Rcpp::wrap(res_nbr_radius),
        Rcpp::Named("pvalue_adj") = Rcpp::wrap(res_pvalue_adj),
        Rcpp::Named("num_nbrs") = Rcpp::wrap(res_num_nbrs),
        Rcpp::Named("expected_num_nbrs") = Rcpp::wrap(res_expected),
        Rcpp::Named("raw_count") = Rcpp::wrap(res_raw_count),
        Rcpp::Named("clump_type") = Rcpp::wrap(res_clump_type),
        Rcpp::Named("is_clumped") = is_clumped,
        Rcpp::Named("all_raw_pvalues") = all_raw_pvalues
    );
}
