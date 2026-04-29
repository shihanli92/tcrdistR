// tcrdist_core.h — Shared types and helpers for all tcrdistR C++ translation units.
//
// Included by every .cpp file in the package.  Contains:
//   - Standard library and Rcpp includes
//   - build_aa_index(): runtime char->column lookup from a named BSD4 matrix
//   - CDR3Data struct and preprocess_cdr3(): precomputed CDR3 metadata
//   - cdr3_dist_fast(): hot-loop CDR3 distance using BSD4_FLAT
//   - VDistLookup struct: named V-gene distance matrix wrapper

#pragma once

#include "blosum.h"

#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------------
// build_aa_index
// ---------------------------------------------------------------------------
// Builds a 26-element lookup array: idx[c - 'A'] = column index in bsd4
// for the 20 standard amino acids, -1 for all others.
//
// Used when BSD4 is received as a named R NumericMatrix (runtime path).
// For the hot constexpr path use tcrdist::aa_to_index() directly.

inline std::array<int, 26> build_aa_index(const Rcpp::NumericMatrix& bsd4) {
    std::array<int, 26> idx;
    idx.fill(-1);

    Rcpp::CharacterVector col_names = Rcpp::colnames(bsd4);
    for (int j = 0; j < col_names.size(); ++j) {
        std::string name = Rcpp::as<std::string>(col_names[j]);
        if (name.size() == 1u) {
            char c = name[0];
            if (c >= 'A' && c <= 'Z') {
                idx[static_cast<int>(c - 'A')] = j;
            }
        }
    }
    return idx;
}

// ---------------------------------------------------------------------------
// CDR3Data — pre-processed CDR3 sequence for hot-loop distance computation
// ---------------------------------------------------------------------------
// Converts amino acid characters to 0-19 uint8_t indices and pre-computes
// gappos / remainder so the inner loop performs only flat array lookups.

struct CDR3Data {
    std::vector<uint8_t> aa;  // AA indices in [0, 19]
    int len;
};

// preprocess_cdr3 — converts a CDR3 string to CDR3Data.
//
// Uses tcrdist::aa_to_index() (constexpr, no runtime matrix needed).
// Throws Rcpp::exception on invalid amino acid characters.

inline CDR3Data preprocess_cdr3(const std::string& seq) {
    CDR3Data d;
    d.len = static_cast<int>(seq.size());
    d.aa.resize(d.len);

    for (int i = 0; i < d.len; ++i) {
        char c = seq[i];
        int idx = tcrdist::aa_to_index(c);
        if (idx < 0) {
            Rcpp::stop("preprocess_cdr3: invalid/unknown amino acid character '%c' in sequence: %s",
                       c, seq.c_str());
        }
        d.aa[i] = static_cast<uint8_t>(idx);
    }

    return d;
}

// ---------------------------------------------------------------------------
// cdr3_dist_fast — core CDR3 distance using BSD4_FLAT constexpr table
// ---------------------------------------------------------------------------
// Matches tcrdist3's nb_tcrdist() with fixed_gappos=False for CDR3:
//   - Same-length: score positions [ntrim, len-ctrim) directly
//   - Different-length: try gap positions from min_gappos to max_gappos,
//     pick the alignment that minimises the substitution score
//
// ntrim = 3, ctrim = 2
// return weight * min_subst_dist + lendiff * gap_penalty

inline double cdr3_dist_fast(const CDR3Data& s1,
                               const CDR3Data& s2,
                               int weight_cdr3_region,
                               int gap_penalty_cdr3_region) {
    // Orient: shortd has <= residues
    const CDR3Data& shortd = (s1.len <= s2.len) ? s1 : s2;
    const CDR3Data& longd  = (s1.len <= s2.len) ? s2 : s1;

    const int lendiff = longd.len - shortd.len;
    const double gap_cost = static_cast<double>(lendiff * gap_penalty_cdr3_region);

    // Sentinel: CDR3 too short to align
    if (shortd.len < 5) {
        return static_cast<double>(weight_cdr3_region * 4.0 * 5 + gap_cost);
    }

    static const int ntrim = 3;
    static const int ctrim = 2;

    const uint8_t* short_aa  = shortd.aa.data();
    const uint8_t* long_aa   = longd.aa.data();
    const int      short_last = shortd.len - 1;
    const int      long_last  = longd.len  - 1;

    // Same-length: no gap needed
    if (lendiff == 0) {
        double dist = 0.0;
        for (int i = ntrim; i < shortd.len - ctrim; ++i) {
            dist += tcrdist::BSD4_FLAT[short_aa[i] * 20 + long_aa[i]];
        }
        return weight_cdr3_region * dist;
    }

    // Different-length: try multiple gap positions, pick minimum
    int min_gappos = 5;
    int max_gappos = shortd.len - 1 - 4;
    while (min_gappos > max_gappos) {
        --min_gappos;
        ++max_gappos;
    }

    double best_dist = -1.0;
    for (int gappos = min_gappos; gappos <= max_gappos; ++gappos) {
        double tmp_dist = 0.0;
        const int remainder = shortd.len - gappos;

        // N-terminal flank: positions [ntrim, gappos)
        for (int i = ntrim; i < gappos; ++i) {
            tmp_dist += tcrdist::BSD4_FLAT[short_aa[i] * 20 + long_aa[i]];
        }
        // C-terminal flank: positions [ctrim, remainder)
        for (int i = ctrim; i < remainder; ++i) {
            tmp_dist += tcrdist::BSD4_FLAT[short_aa[short_last - i] * 20
                                         + long_aa[long_last  - i]];
        }

        if (tmp_dist < best_dist || best_dist < 0.0) {
            best_dist = tmp_dist;
        }
        if (best_dist == 0.0) break;
    }

    return weight_cdr3_region * best_dist + gap_cost;
}

// ---------------------------------------------------------------------------
// cdr3_dist_fast (budgeted overload) — early termination for background
// ---------------------------------------------------------------------------
// Same algorithm, but returns early when gap_cost alone exceeds `budget`.
// Used in background distribution computation where v_dist is already known
// and the caller passes budget = max_dist - v_dist.

inline double cdr3_dist_fast(const CDR3Data& s1,
                               const CDR3Data& s2,
                               int weight_cdr3_region,
                               int gap_penalty_cdr3_region,
                               double budget) {
    const CDR3Data& shortd = (s1.len <= s2.len) ? s1 : s2;
    const CDR3Data& longd  = (s1.len <= s2.len) ? s2 : s1;

    const int lendiff = longd.len - shortd.len;
    const double gap_cost = static_cast<double>(lendiff * gap_penalty_cdr3_region);

    if (shortd.len < 5) {
        return static_cast<double>(weight_cdr3_region * 4.0 * 5 + gap_cost);
    }

    // Early exit: gap penalty alone exceeds remaining budget
    if (gap_cost >= budget) {
        return gap_cost;
    }

    static const int ntrim = 3;
    static const int ctrim = 2;

    const uint8_t* short_aa  = shortd.aa.data();
    const uint8_t* long_aa   = longd.aa.data();
    const int      short_last = shortd.len - 1;
    const int      long_last  = longd.len  - 1;

    // Remaining budget after gap cost, expressed as raw substitution limit
    const double subst_budget = (budget - gap_cost) / weight_cdr3_region;

    if (lendiff == 0) {
        double dist = 0.0;
        for (int i = ntrim; i < shortd.len - ctrim; ++i) {
            dist += tcrdist::BSD4_FLAT[short_aa[i] * 20 + long_aa[i]];
            if (dist > subst_budget) return weight_cdr3_region * dist;
        }
        return weight_cdr3_region * dist;
    }

    int min_gappos = 5;
    int max_gappos = shortd.len - 1 - 4;
    while (min_gappos > max_gappos) {
        --min_gappos;
        ++max_gappos;
    }

    double best_dist = -1.0;
    for (int gappos = min_gappos; gappos <= max_gappos; ++gappos) {
        double tmp_dist = 0.0;
        const int remainder = shortd.len - gappos;

        for (int i = ntrim; i < gappos; ++i) {
            tmp_dist += tcrdist::BSD4_FLAT[short_aa[i] * 20 + long_aa[i]];
        }
        for (int i = ctrim; i < remainder; ++i) {
            tmp_dist += tcrdist::BSD4_FLAT[short_aa[short_last - i] * 20
                                         + long_aa[long_last  - i]];
        }

        if (tmp_dist < best_dist || best_dist < 0.0) {
            best_dist = tmp_dist;
        }
        if (best_dist == 0.0) break;
    }

    return weight_cdr3_region * best_dist + gap_cost;
}

// ---------------------------------------------------------------------------
// VDistLookup — named V-gene distance matrix wrapper
// ---------------------------------------------------------------------------
// Builds a gene-name -> index map from the row/column names of a named
// NumericMatrix, then stores a flat row-major copy of the data for O(1)
// lookups in the inner loop without going through Rcpp overhead.

struct VDistLookup {
    std::unordered_map<std::string, int> gene_to_idx;
    std::vector<double>                  dist_flat;   // row-major copy
    int                                  n_genes = 0;

    // Build from a named symmetric NumericMatrix (row names = gene IDs).
    void build(const Rcpp::NumericMatrix& v_dist_mat) {
        Rcpp::CharacterVector rnames = Rcpp::rownames(v_dist_mat);
        n_genes = static_cast<int>(rnames.size());
        gene_to_idx.reserve(n_genes);
        for (int k = 0; k < n_genes; ++k) {
            gene_to_idx[Rcpp::as<std::string>(rnames[k])] = k;
        }

        // Copy column-major R matrix to row-major C++ vector.
        // R element [r, c] is stored at c * nrow + r (column-major).
        // We want dist_flat[r * n_genes + c] = mat[r, c].
        dist_flat.resize(static_cast<size_t>(n_genes) * n_genes);
        const double* ptr = v_dist_mat.begin();
        for (int c = 0; c < n_genes; ++c) {
            for (int r = 0; r < n_genes; ++r) {
                dist_flat[r * n_genes + c] = ptr[c * n_genes + r];
            }
        }
    }

    // Lookup by pre-resolved row/column indices — O(1).
    inline double lookup(int i, int j) const {
        return dist_flat[i * n_genes + j];
    }

    // Lookup by gene name strings — O(1) average (hash map).
    double lookup(const std::string& g1, const std::string& g2) const {
        auto it1 = gene_to_idx.find(g1);
        auto it2 = gene_to_idx.find(g2);
        if (it1 == gene_to_idx.end()) {
            Rcpp::stop("VDistLookup::lookup: gene '%s' not found in V-region matrix",
                       g1.c_str());
        }
        if (it2 == gene_to_idx.end()) {
            Rcpp::stop("VDistLookup::lookup: gene '%s' not found in V-region matrix",
                       g2.c_str());
        }
        return lookup(it1->second, it2->second);
    }

    // Resolve a gene name to its row/column index, returning -1 if not found.
    int resolve(const std::string& gene) const {
        auto it = gene_to_idx.find(gene);
        return (it != gene_to_idx.end()) ? it->second : -1;
    }
};

// ---------------------------------------------------------------------------
// PreparedTCRs — pre-processed TCR set for distance computation
// ---------------------------------------------------------------------------
// Converts R CharacterVectors to C++ types, resolves V-gene names to integer
// indices, and preprocesses CDR3 sequences — the boilerplate shared by
// tcrdist_matrix, tcrdist_sparse, tcrdist_knn, tcrdist_radius, tcrdist_rect.

struct PreparedTCRs {
    int                   n;
    std::vector<int>      vi_a;   // resolved V-gene alpha indices
    std::vector<int>      vi_b;   // resolved V-gene beta indices
    std::vector<CDR3Data> cdr3a;  // preprocessed CDR3 alpha
    std::vector<CDR3Data> cdr3b;  // preprocessed CDR3 beta
};

// prepare_tcrs — converts R inputs to PreparedTCRs.
//
// Validates that all four CharacterVectors have the same length, converts
// strings, resolves V-gene names to integer indices (fails fast on unknown
// genes), and preprocesses CDR3 sequences.  `caller` is used in error messages.

inline PreparedTCRs prepare_tcrs(
    const Rcpp::CharacterVector& va_genes,
    const Rcpp::CharacterVector& cdr3a_seqs,
    const Rcpp::CharacterVector& vb_genes,
    const Rcpp::CharacterVector& cdr3b_seqs,
    const VDistLookup& vla,
    const VDistLookup& vlb,
    const char* caller
) {
    const int n = va_genes.size();

    if (cdr3a_seqs.size() != n || vb_genes.size() != n || cdr3b_seqs.size() != n) {
        Rcpp::stop(
            "%s: va_genes, cdr3a_seqs, vb_genes, cdr3b_seqs "
            "must all have the same length", caller
        );
    }

    PreparedTCRs t;
    t.n = n;
    t.vi_a.resize(n);
    t.vi_b.resize(n);
    t.cdr3a.resize(n);
    t.cdr3b.resize(n);

    for (int k = 0; k < n; ++k) {
        std::string va_str  = Rcpp::as<std::string>(va_genes[k]);
        std::string vb_str  = Rcpp::as<std::string>(vb_genes[k]);

        int ia = vla.resolve(va_str);
        if (ia < 0) {
            Rcpp::stop("%s: alpha V-gene '%s' not found in v_dist_a",
                       caller, va_str.c_str());
        }
        t.vi_a[k] = ia;

        int ib = vlb.resolve(vb_str);
        if (ib < 0) {
            Rcpp::stop("%s: beta V-gene '%s' not found in v_dist_b",
                       caller, vb_str.c_str());
        }
        t.vi_b[k] = ib;

        t.cdr3a[k] = preprocess_cdr3(Rcpp::as<std::string>(cdr3a_seqs[k]));
        t.cdr3b[k] = preprocess_cdr3(Rcpp::as<std::string>(cdr3b_seqs[k]));
    }

    return t;
}
