// tcrdist_background.cpp — Background distribution computation for TCR clumping
//
// For each foreground paired-chain TCR, computes the cumulative frequency
// distribution of TCRdist distances to all background chain pairs.  This is
// the inner engine of .estimate_background_tcrdist_distributions() moved
// entirely to C++ to eliminate millions of R-to-C++ call transitions.
//
// Extracted from rconga/src/tcrdist_rcpp.cpp lines 630-942, refactored to
// use tcrdistR shared infrastructure (CDR3Data, preprocess_cdr3, VDistLookup,
// constexpr BSD4_FLAT, budgeted cdr3_dist_fast overload).

#include "tcrdist_core.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_calc_background_distributions(
    const Rcpp::CharacterVector& fg_va,
    const Rcpp::CharacterVector& fg_cdr3a,
    const Rcpp::CharacterVector& fg_vb,
    const Rcpp::CharacterVector& fg_cdr3b,
    const Rcpp::CharacterVector& bg_va,
    const Rcpp::CharacterVector& bg_cdr3a,
    const Rcpp::CharacterVector& bg_vb,
    const Rcpp::CharacterVector& bg_cdr3b,
    const Rcpp::NumericMatrix& v_dist_a,
    const Rcpp::NumericMatrix& v_dist_b,
    int max_dist,
    double pseudocount = 0.25,
    int weight_cdr3_region = 3,
    int gap_penalty_cdr3_region = 12
) {
    const int n_fg   = fg_va.size();
    const int n_bg_a = bg_va.size();
    const int n_bg_b = bg_vb.size();
    const int n_bins = max_dist + 1;
    const double n_bg_pairs = static_cast<double>(n_bg_a) * n_bg_b;

    // ---- validate input lengths ---------------------------------------------
    if (fg_cdr3a.size() != n_fg || fg_vb.size() != n_fg || fg_cdr3b.size() != n_fg) {
        Rcpp::stop(
            "rcpp_calc_background_distributions: fg_va, fg_cdr3a, fg_vb, fg_cdr3b "
            "must all have the same length (%d vs %d/%d/%d)",
            n_fg, static_cast<int>(fg_cdr3a.size()),
            static_cast<int>(fg_vb.size()), static_cast<int>(fg_cdr3b.size()));
    }
    if (bg_cdr3a.size() != n_bg_a) {
        Rcpp::stop(
            "rcpp_calc_background_distributions: bg_va and bg_cdr3a must have "
            "the same length (%d vs %d)", n_bg_a, static_cast<int>(bg_cdr3a.size()));
    }
    if (bg_cdr3b.size() != n_bg_b) {
        Rcpp::stop(
            "rcpp_calc_background_distributions: bg_vb and bg_cdr3b must have "
            "the same length (%d vs %d)", n_bg_b, static_cast<int>(bg_cdr3b.size()));
    }

    // ---- build V-gene distance lookups using VDistLookup --------------------
    VDistLookup vlookup_a, vlookup_b;
    vlookup_a.build(v_dist_a);
    vlookup_b.build(v_dist_b);

    // ---- pre-convert CDR3 strings to indexed form ---------------------------
    // preprocess_cdr3() from tcrdist_core.h uses constexpr aa_to_index()
    // and throws on invalid characters.  We warn for short sequences since
    // cdr3_dist_fast() returns a sentinel distance for them.

    std::vector<CDR3Data> fg_cdr3a_d(n_fg), fg_cdr3b_d(n_fg);
    for (int k = 0; k < n_fg; ++k) {
        std::string seq_a = Rcpp::as<std::string>(fg_cdr3a[k]);
        std::string seq_b = Rcpp::as<std::string>(fg_cdr3b[k]);
        if (static_cast<int>(seq_a.size()) < 5)
            Rf_warning("rcpp_calc_background_distributions: fg CDR3a too short (got %d): '%s'",
                       static_cast<int>(seq_a.size()), seq_a.c_str());
        if (static_cast<int>(seq_b.size()) < 5)
            Rf_warning("rcpp_calc_background_distributions: fg CDR3b too short (got %d): '%s'",
                       static_cast<int>(seq_b.size()), seq_b.c_str());
        fg_cdr3a_d[k] = preprocess_cdr3(seq_a);
        fg_cdr3b_d[k] = preprocess_cdr3(seq_b);
    }

    std::vector<CDR3Data> bg_cdr3a_d(n_bg_a), bg_cdr3b_d(n_bg_b);
    for (int k = 0; k < n_bg_a; ++k) {
        std::string seq = Rcpp::as<std::string>(bg_cdr3a[k]);
        if (static_cast<int>(seq.size()) < 5)
            Rf_warning("rcpp_calc_background_distributions: bg CDR3a too short (got %d): '%s'",
                       static_cast<int>(seq.size()), seq.c_str());
        bg_cdr3a_d[k] = preprocess_cdr3(seq);
    }
    for (int k = 0; k < n_bg_b; ++k) {
        std::string seq = Rcpp::as<std::string>(bg_cdr3b[k]);
        if (static_cast<int>(seq.size()) < 5)
            Rf_warning("rcpp_calc_background_distributions: bg CDR3b too short (got %d): '%s'",
                       static_cast<int>(seq.size()), seq.c_str());
        bg_cdr3b_d[k] = preprocess_cdr3(seq);
    }

    // ---- pre-resolve V-gene indices for foreground chains -------------------
    std::vector<int> fg_a_vidx(n_fg), fg_b_vidx(n_fg);
    for (int k = 0; k < n_fg; ++k) {
        fg_a_vidx[k] = vlookup_a.resolve(Rcpp::as<std::string>(fg_va[k]));
        fg_b_vidx[k] = vlookup_b.resolve(Rcpp::as<std::string>(fg_vb[k]));
    }

    // ---- pre-resolve V-gene indices for background chains -------------------
    // -1 means gene not found; V-region distance treated as 0 (matches R).
    std::vector<int> bg_a_vidx(n_bg_a), bg_b_vidx(n_bg_b);
    for (int j = 0; j < n_bg_a; ++j) {
        bg_a_vidx[j] = vlookup_a.resolve(Rcpp::as<std::string>(bg_va[j]));
    }
    for (int j = 0; j < n_bg_b; ++j) {
        bg_b_vidx[j] = vlookup_b.resolve(Rcpp::as<std::string>(bg_vb[j]));
    }

    // ---- allocate result and temp buffers -----------------------------------
    Rcpp::NumericMatrix result(n_fg, n_bins);
    std::vector<double> a_hist(n_bins);
    std::vector<double> b_hist(n_bins);
    std::vector<double> paired(n_bins);

    // ---- main loop: for each foreground TCR ---------------------------------
    for (int ii = 0; ii < n_fg; ++ii) {
        if (ii % 50 == 0) {
            Rcpp::checkUserInterrupt();
        }

        const int fg_a_vi = fg_a_vidx[ii];
        const int fg_b_vi = fg_b_vidx[ii];

        // ---- Alpha distances: fg[ii] vs all bg_alpha ------------------------
        std::fill(a_hist.begin(), a_hist.end(), 0.0);
        const CDR3Data& fg_a = fg_cdr3a_d[ii];

        for (int j = 0; j < n_bg_a; ++j) {
            double v_dist = 0.0;
            if (fg_a_vi >= 0 && bg_a_vidx[j] >= 0) {
                v_dist = vlookup_a.lookup(fg_a_vi, bg_a_vidx[j]);
            }
            // Early termination: V-distance alone >= max_dist
            if (v_dist >= max_dist) {
                a_hist[max_dist] += 1.0;
                continue;
            }
            double d = v_dist + cdr3_dist_fast(fg_a, bg_cdr3a_d[j],
                                               weight_cdr3_region,
                                               gap_penalty_cdr3_region,
                                               max_dist - v_dist);
            int bin = static_cast<int>(d);
            if (bin > max_dist) bin = max_dist;
            if (bin < 0) bin = 0;
            a_hist[bin] += 1.0;
        }

        // ---- Beta distances: fg[ii] vs all bg_beta -------------------------
        std::fill(b_hist.begin(), b_hist.end(), 0.0);
        const CDR3Data& fg_b = fg_cdr3b_d[ii];

        for (int j = 0; j < n_bg_b; ++j) {
            double v_dist = 0.0;
            if (fg_b_vi >= 0 && bg_b_vidx[j] >= 0) {
                v_dist = vlookup_b.lookup(fg_b_vi, bg_b_vidx[j]);
            }
            if (v_dist >= max_dist) {
                b_hist[max_dist] += 1.0;
                continue;
            }
            double d = v_dist + cdr3_dist_fast(fg_b, bg_cdr3b_d[j],
                                               weight_cdr3_region,
                                               gap_penalty_cdr3_region,
                                               max_dist - v_dist);
            int bin = static_cast<int>(d);
            if (bin > max_dist) bin = max_dist;
            if (bin < 0) bin = 0;
            b_hist[bin] += 1.0;
        }

        // ---- Convolve: paired[d] = sum_{k=0}^{d} a[k] * b[d-k] ------------
        std::fill(paired.begin(), paired.end(), 0.0);
        for (int d = 0; d < n_bins; ++d) {
            double sum = 0.0;
            for (int k = 0; k <= d; ++k) {
                sum += a_hist[k] * b_hist[d - k];
            }
            paired[d] = sum;
        }

        // ---- Cumulative sum + pseudocount + normalize -----------------------
        double cum = 0.0;
        for (int d = 0; d < n_bins; ++d) {
            cum += paired[d];
            double val = (cum > pseudocount) ? cum : pseudocount;
            result(ii, d) = val / n_bg_pairs;
        }
    }

    return result;
}
