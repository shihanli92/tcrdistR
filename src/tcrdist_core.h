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
    int gappos;     // gap insertion point (0 for short seqs)
    int remainder;  // len - gappos (0 for short seqs)
};

// preprocess_cdr3 — converts a CDR3 string to CDR3Data.
//
// Uses tcrdist::aa_to_index() (constexpr, no runtime matrix needed).
// For sequences with len < 5, gappos and remainder are set to 0; the
// sentinel return in cdr3_dist_fast() handles them correctly.
//
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

    if (d.len >= 5) {
        d.gappos    = std::min(6, 3 + (d.len - 5) / 2);
        d.remainder = d.len - d.gappos;
    } else {
        d.gappos    = 0;
        d.remainder = 0;
    }
    return d;
}

// ---------------------------------------------------------------------------
// cdr3_dist_fast — core CDR3 distance using BSD4_FLAT constexpr table
// ---------------------------------------------------------------------------
// Direct port of rconga/src/tcrdist_rcpp.cpp lines 702-741 (cdr3_dist_fast).
// Uses the constexpr BSD4_FLAT array rather than a passed-in matrix pointer,
// so no setup is required at call sites.
//
// Algorithm (ALIGN_CDR3S = FALSE, TRIM_CDR3S = TRUE):
//   ntrim = 3, ctrim = 2
//   gappos = min(6, 3 + (lenshort - 5) / 2)   [integer division]
//   remainder = lenshort - gappos
//   N-terminal: sum BSD4[short[i]][long[i]] for i in [ntrim, gappos)
//   C-terminal: sum BSD4[short[lenshort-1-i]][long[lenlong-1-i]] for i in [ctrim, remainder)
//   return weight * dist + lendiff * gap_penalty

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

    double dist = 0.0;

    // N-terminal flank: positions [ntrim, gappos)
    for (int i = ntrim; i < shortd.gappos; ++i) {
        dist += tcrdist::BSD4_FLAT[short_aa[i] * 20 + long_aa[i]];
    }

    // C-terminal flank: positions [ctrim, remainder)
    for (int i = ctrim; i < shortd.remainder; ++i) {
        dist += tcrdist::BSD4_FLAT[short_aa[short_last - i] * 20
                                 + long_aa[long_last  - i]];
    }

    return weight_cdr3_region * dist + gap_cost;
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

    double dist = 0.0;

    for (int i = ntrim; i < shortd.gappos; ++i) {
        dist += tcrdist::BSD4_FLAT[short_aa[i] * 20 + long_aa[i]];
    }

    for (int i = ctrim; i < shortd.remainder; ++i) {
        dist += tcrdist::BSD4_FLAT[short_aa[short_last - i] * 20
                                 + long_aa[long_last  - i]];
    }

    return weight_cdr3_region * dist + gap_cost;
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
