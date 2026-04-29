// tcr_sampler.cpp — TCR sampling and allele-optimization helpers.
//
// Exports:
//   rcpp_count_nuc_matches()               — prefix nucleotide match scoring
//   rcpp_find_alternate_alleles_batch()     — batch allele optimization
//   rcpp_resample_shuffled_tcr_chains()     — junction breakpoint splicing
//
// Direct port from rconga/src/tcr_sampler_rcpp.cpp (512 lines).
// Algorithms are identical; only packaging changes:
//   - includes tcrdist_core.h (package shared header)
//   - inline codon table replaces passed-in codon_keys/codon_values
//   - error messages use function-name prefixes

#include "tcrdist_core.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace Rcpp;


// ---------------------------------------------------------------------------
// 1. rcpp_count_nuc_matches — Direct port of .count_nuc_matches()
//
// Called 500K+ times from .find_alternate_alleles() and .analyze_junction().
// The R version does strsplit() + character-by-character loop; this version
// uses direct std::string indexing with zero allocation overhead.
// ---------------------------------------------------------------------------

//' Count nucleotide prefix matches between two sequences (C++ implementation)
//'
//' Linear scan prefix matching: +1 per match, mismatch_score per mismatch.
//' Returns the length of the best-scoring prefix (the position at which the
//' cumulative score was highest).
//'
//' @param a Character string. First nucleotide sequence.
//' @param b Character string. Second nucleotide sequence.
//' @param mismatch_score Integer. Score added per mismatch position.
//'   Default -4.
//' @return Integer: length of best-scoring prefix (0 if either string is
//'   empty).
//' @examples
//' \dontrun{
//'   rcpp_count_nuc_matches("ATCGATCG", "ATCGTTCG")
//' }
//' @export
// [[Rcpp::export]]
int rcpp_count_nuc_matches(const std::string& a, const std::string& b,
                           int mismatch_score = -4) {
    int na = static_cast<int>(a.size());
    int nb = static_cast<int>(b.size());
    if (na == 0 || nb == 0) return 0;

    int n = std::min(na, nb);
    int match_score = 1;
    int best_score  = 0;
    int score       = 0;
    int num_matches = 0;

    for (int i = 0; i < n; ++i) {
        if (a[i] == b[i]) {
            score += match_score;
        } else {
            score += mismatch_score;
        }
        if (score >= best_score) {
            best_score  = score;
            num_matches = i + 1;  // 1-based count
        }
    }

    return num_matches;
}


// ---------------------------------------------------------------------------
// 2. rcpp_find_alternate_alleles_batch — Batch allele optimization
//
// Moves the entire .find_alternate_alleles_for_tcrs() inner loop into C++.
// For each TCR, iterates over V and J alternate alleles, performing
// count_nuc_matches (with string reversal for J alleles) entirely in C++.
// ---------------------------------------------------------------------------

namespace {

// Inline count_nuc_matches for use within this translation unit
int count_nuc_matches_inline(const std::string& a, const std::string& b,
                             int mismatch_score) {
    int n = std::min(static_cast<int>(a.size()), static_cast<int>(b.size()));
    if (n == 0) return 0;
    int best_score = 0, score = 0, num_matches = 0;
    for (int i = 0; i < n; ++i) {
        score += (a[i] == b[i]) ? 1 : mismatch_score;
        if (score >= best_score) {
            best_score = score;
            num_matches = i + 1;
        }
    }
    return num_matches;
}

}  // anonymous namespace


//' Find alternate alleles for a batch of TCRs (C++ implementation)
//'
//' For each paired TCR (alpha+beta chains), tries alternate V/J alleles
//' (same gene, different allele number after \code{*}) and picks the allele
//' whose CDR3 nucleotide prefix match score improves by at least
//' \code{min_improvement} over the current allele.  J alleles are compared
//' using reversed sequences.
//'
//' @param va_genes Character vector. Alpha V gene names.
//' @param ja_genes Character vector. Alpha J gene names.
//' @param cdr3a_nucseqs Character vector. Alpha CDR3 nucleotide sequences.
//' @param vb_genes Character vector. Beta V gene names.
//' @param jb_genes Character vector. Beta J gene names.
//' @param cdr3b_nucseqs Character vector. Beta CDR3 nucleotide sequences.
//' @param all_gene_names Character vector. All gene names in the database.
//' @param all_v_cdr3_nucseqs Character vector. Pre-computed V CDR3 nucleotide
//'   sequences (same order as \code{all_gene_names}).
//' @param all_j_cdr3_nucseqs Character vector. Pre-computed J CDR3 nucleotide
//'   sequences (same order as \code{all_gene_names}).
//' @param mismatch_score Integer. Mismatch penalty for prefix scoring.
//'   Default -4.
//' @param min_improvement Integer. Minimum score improvement required to
//'   accept an alternate allele over the current one. Default 2.
//' @return A list with elements: \code{new_va}, \code{new_ja}, \code{new_vb},
//'   \code{new_jb} (character vectors), and \code{counts} (data.frame with
//'   columns \code{old_gene}, \code{new_gene}, \code{count}).
//' @examples
//' \dontrun{
//'   result <- rcpp_find_alternate_alleles_batch(
//'     va_genes, ja_genes, cdr3a_nucseqs,
//'     vb_genes, jb_genes, cdr3b_nucseqs,
//'     all_gene_names, all_v_cdr3_nucseqs, all_j_cdr3_nucseqs
//'   )
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_find_alternate_alleles_batch(
    const CharacterVector& va_genes,
    const CharacterVector& ja_genes,
    const CharacterVector& cdr3a_nucseqs,
    const CharacterVector& vb_genes,
    const CharacterVector& jb_genes,
    const CharacterVector& cdr3b_nucseqs,
    const CharacterVector& all_gene_names,
    const CharacterVector& all_v_cdr3_nucseqs,
    const CharacterVector& all_j_cdr3_nucseqs,
    int mismatch_score = -4,
    int min_improvement = 2) {

    int n = va_genes.size();

    // Pre-convert R CharacterVectors to C++ strings for fast access
    int n_genes = all_gene_names.size();
    std::vector<std::string> gene_names(n_genes);
    std::vector<std::string> v_nucseqs(n_genes);
    std::vector<std::string> j_nucseqs(n_genes);

    for (int i = 0; i < n_genes; ++i) {
        gene_names[i] = Rcpp::as<std::string>(all_gene_names[i]);
        v_nucseqs[i]  = Rcpp::as<std::string>(all_v_cdr3_nucseqs[i]);
        j_nucseqs[i]  = Rcpp::as<std::string>(all_j_cdr3_nucseqs[i]);
    }

    // Build prefix -> gene indices map for alternate allele lookup
    // prefix = gene name up to and including "*"
    std::unordered_map<std::string, std::vector<int>> prefix_to_indices;
    for (int i = 0; i < n_genes; ++i) {
        size_t star_pos = gene_names[i].find('*');
        if (star_pos != std::string::npos) {
            std::string prefix = gene_names[i].substr(0, star_pos + 1);
            prefix_to_indices[prefix].push_back(i);
        }
    }

    // Build name -> index map
    std::unordered_map<std::string, int> name_to_idx;
    for (int i = 0; i < n_genes; ++i) {
        name_to_idx[gene_names[i]] = i;
    }

    // Output vectors
    CharacterVector new_va(n), new_ja(n), new_vb(n), new_jb(n);

    // Also track counts for the allowed-swaps logic
    // counts[old_gene][new_gene] = count
    std::unordered_map<std::string, std::unordered_map<std::string, int>> all_counts;

    for (int ii = 0; ii < n; ++ii) {
        std::string cur_va = Rcpp::as<std::string>(va_genes[ii]);
        std::string cur_ja = Rcpp::as<std::string>(ja_genes[ii]);
        std::string cdr3a  = Rcpp::as<std::string>(cdr3a_nucseqs[ii]);
        std::string cur_vb = Rcpp::as<std::string>(vb_genes[ii]);
        std::string cur_jb = Rcpp::as<std::string>(jb_genes[ii]);
        std::string cdr3b  = Rcpp::as<std::string>(cdr3b_nucseqs[ii]);

        // --- Alpha V alleles ---
        std::string best_va = cur_va;
        int best_va_matched = 0;
        {
            size_t star_pos = cur_va.find('*');
            std::string prefix = (star_pos != std::string::npos) ?
                cur_va.substr(0, star_pos + 1) : cur_va;
            auto it = prefix_to_indices.find(prefix);
            if (it != prefix_to_indices.end()) {
                // First check the current gene
                auto idx_it = name_to_idx.find(cur_va);
                if (idx_it != name_to_idx.end() &&
                    !v_nucseqs[idx_it->second].empty()) {
                    best_va_matched = count_nuc_matches_inline(
                        v_nucseqs[idx_it->second], cdr3a, mismatch_score);
                }
                // Then check alternates
                for (int gi : it->second) {
                    if (gene_names[gi] == cur_va) continue;
                    if (v_nucseqs[gi].empty()) continue;
                    int matched = count_nuc_matches_inline(
                        v_nucseqs[gi], cdr3a, mismatch_score);
                    if (matched > best_va_matched &&
                        (best_va == cur_va ?
                         matched - best_va_matched >= min_improvement : true)) {
                        best_va = gene_names[gi];
                        best_va_matched = matched;
                    }
                }
            }
        }

        // --- Alpha J alleles ---
        std::string best_ja = cur_ja;
        int best_ja_matched = 0;
        {
            std::string cdr3a_rev(cdr3a.rbegin(), cdr3a.rend());
            size_t star_pos = cur_ja.find('*');
            std::string prefix = (star_pos != std::string::npos) ?
                cur_ja.substr(0, star_pos + 1) : cur_ja;
            auto it = prefix_to_indices.find(prefix);
            if (it != prefix_to_indices.end()) {
                auto idx_it = name_to_idx.find(cur_ja);
                if (idx_it != name_to_idx.end() &&
                    !j_nucseqs[idx_it->second].empty()) {
                    std::string j_rev(j_nucseqs[idx_it->second].rbegin(),
                                      j_nucseqs[idx_it->second].rend());
                    best_ja_matched = count_nuc_matches_inline(
                        j_rev, cdr3a_rev, mismatch_score);
                }
                for (int gi : it->second) {
                    if (gene_names[gi] == cur_ja) continue;
                    if (j_nucseqs[gi].empty()) continue;
                    std::string j_rev(j_nucseqs[gi].rbegin(),
                                      j_nucseqs[gi].rend());
                    int matched = count_nuc_matches_inline(
                        j_rev, cdr3a_rev, mismatch_score);
                    if (matched > best_ja_matched &&
                        (best_ja == cur_ja ?
                         matched - best_ja_matched >= min_improvement : true)) {
                        best_ja = gene_names[gi];
                        best_ja_matched = matched;
                    }
                }
            }
        }

        // --- Beta V alleles ---
        std::string best_vb = cur_vb;
        int best_vb_matched = 0;
        {
            size_t star_pos = cur_vb.find('*');
            std::string prefix = (star_pos != std::string::npos) ?
                cur_vb.substr(0, star_pos + 1) : cur_vb;
            auto it = prefix_to_indices.find(prefix);
            if (it != prefix_to_indices.end()) {
                auto idx_it = name_to_idx.find(cur_vb);
                if (idx_it != name_to_idx.end() &&
                    !v_nucseqs[idx_it->second].empty()) {
                    best_vb_matched = count_nuc_matches_inline(
                        v_nucseqs[idx_it->second], cdr3b, mismatch_score);
                }
                for (int gi : it->second) {
                    if (gene_names[gi] == cur_vb) continue;
                    if (v_nucseqs[gi].empty()) continue;
                    int matched = count_nuc_matches_inline(
                        v_nucseqs[gi], cdr3b, mismatch_score);
                    if (matched > best_vb_matched &&
                        (best_vb == cur_vb ?
                         matched - best_vb_matched >= min_improvement : true)) {
                        best_vb = gene_names[gi];
                        best_vb_matched = matched;
                    }
                }
            }
        }

        // --- Beta J alleles ---
        std::string best_jb = cur_jb;
        int best_jb_matched = 0;
        {
            std::string cdr3b_rev(cdr3b.rbegin(), cdr3b.rend());
            size_t star_pos = cur_jb.find('*');
            std::string prefix = (star_pos != std::string::npos) ?
                cur_jb.substr(0, star_pos + 1) : cur_jb;
            auto it = prefix_to_indices.find(prefix);
            if (it != prefix_to_indices.end()) {
                auto idx_it = name_to_idx.find(cur_jb);
                if (idx_it != name_to_idx.end() &&
                    !j_nucseqs[idx_it->second].empty()) {
                    std::string j_rev(j_nucseqs[idx_it->second].rbegin(),
                                      j_nucseqs[idx_it->second].rend());
                    best_jb_matched = count_nuc_matches_inline(
                        j_rev, cdr3b_rev, mismatch_score);
                }
                for (int gi : it->second) {
                    if (gene_names[gi] == cur_jb) continue;
                    if (j_nucseqs[gi].empty()) continue;
                    std::string j_rev(j_nucseqs[gi].rbegin(),
                                      j_nucseqs[gi].rend());
                    int matched = count_nuc_matches_inline(
                        j_rev, cdr3b_rev, mismatch_score);
                    if (matched > best_jb_matched &&
                        (best_jb == cur_jb ?
                         matched - best_jb_matched >= min_improvement : true)) {
                        best_jb = gene_names[gi];
                        best_jb_matched = matched;
                    }
                }
            }
        }

        new_va[ii] = best_va;
        new_ja[ii] = best_ja;
        new_vb[ii] = best_vb;
        new_jb[ii] = best_jb;

        // Track counts
        all_counts[cur_va][best_va]++;
        all_counts[cur_ja][best_ja]++;
        all_counts[cur_vb][best_vb]++;
        all_counts[cur_jb][best_jb]++;
    }

    // Flatten all_counts into a DataFrame (old_gene, new_gene, count)
    std::vector<std::string> cnt_old, cnt_new;
    std::vector<int> cnt_val;
    for (const auto& outer : all_counts) {
        for (const auto& inner : outer.second) {
            cnt_old.push_back(outer.first);
            cnt_new.push_back(inner.first);
            cnt_val.push_back(inner.second);
        }
    }

    Rcpp::DataFrame counts_df = Rcpp::DataFrame::create(
        Rcpp::Named("old_gene") = Rcpp::wrap(cnt_old),
        Rcpp::Named("new_gene") = Rcpp::wrap(cnt_new),
        Rcpp::Named("count")    = Rcpp::wrap(cnt_val),
        Rcpp::Named("stringsAsFactors") = false
    );

    return Rcpp::List::create(
        Rcpp::Named("new_va") = new_va,
        Rcpp::Named("new_ja") = new_ja,
        Rcpp::Named("new_vb") = new_vb,
        Rcpp::Named("new_jb") = new_jb,
        Rcpp::Named("counts") = counts_df
    );
}


// ---------------------------------------------------------------------------
// 3. rcpp_resample_shuffled_tcr_chains — Junction breakpoint splicing
//
// Moves the while-loop from .resample_shuffled_tcr_chains() into C++.
// Key speedups: std::string::substr() replaces substring()+paste0(),
// inline codon translation replaces get_translation(), and stop codon
// detection is done during translation (no grepl).
// ---------------------------------------------------------------------------

namespace {

// Standard genetic code: 64 codons -> amino acid (or '*' for stop).
// Built once per call via init_codon_table().
std::unordered_map<std::string, char> init_codon_table() {
    std::unordered_map<std::string, char> ct;
    ct.reserve(64);
    // Phe
    ct["TTT"] = 'F'; ct["TTC"] = 'F';
    // Leu
    ct["TTA"] = 'L'; ct["TTG"] = 'L';
    ct["CTT"] = 'L'; ct["CTC"] = 'L'; ct["CTA"] = 'L'; ct["CTG"] = 'L';
    // Ile
    ct["ATT"] = 'I'; ct["ATC"] = 'I'; ct["ATA"] = 'I';
    // Met (start)
    ct["ATG"] = 'M';
    // Val
    ct["GTT"] = 'V'; ct["GTC"] = 'V'; ct["GTA"] = 'V'; ct["GTG"] = 'V';
    // Ser
    ct["TCT"] = 'S'; ct["TCC"] = 'S'; ct["TCA"] = 'S'; ct["TCG"] = 'S';
    ct["AGT"] = 'S'; ct["AGC"] = 'S';
    // Pro
    ct["CCT"] = 'P'; ct["CCC"] = 'P'; ct["CCA"] = 'P'; ct["CCG"] = 'P';
    // Thr
    ct["ACT"] = 'T'; ct["ACC"] = 'T'; ct["ACA"] = 'T'; ct["ACG"] = 'T';
    // Ala
    ct["GCT"] = 'A'; ct["GCC"] = 'A'; ct["GCA"] = 'A'; ct["GCG"] = 'A';
    // Tyr
    ct["TAT"] = 'Y'; ct["TAC"] = 'Y';
    // Stop
    ct["TAA"] = '*'; ct["TAG"] = '*'; ct["TGA"] = '*';
    // His
    ct["CAT"] = 'H'; ct["CAC"] = 'H';
    // Gln
    ct["CAA"] = 'Q'; ct["CAG"] = 'Q';
    // Asn
    ct["AAT"] = 'N'; ct["AAC"] = 'N';
    // Lys
    ct["AAA"] = 'K'; ct["AAG"] = 'K';
    // Asp
    ct["GAT"] = 'D'; ct["GAC"] = 'D';
    // Glu
    ct["GAA"] = 'E'; ct["GAG"] = 'E';
    // Cys
    ct["TGT"] = 'C'; ct["TGC"] = 'C';
    // Trp
    ct["TGG"] = 'W';
    // Arg
    ct["CGT"] = 'R'; ct["CGC"] = 'R'; ct["CGA"] = 'R'; ct["CGG"] = 'R';
    ct["AGA"] = 'R'; ct["AGG"] = 'R';
    // Gly
    ct["GGT"] = 'G'; ct["GGC"] = 'G'; ct["GGA"] = 'G'; ct["GGG"] = 'G';

    return ct;
}

}  // anonymous namespace


//' Resample shuffled TCR chains by junction breakpoint splicing (C++ implementation)
//'
//' Randomly selects pairs of junctions, finds shared breakpoints, and splices
//' chimeric nucleotide sequences.  Validates that the result has a length
//' divisible by 3 and contains no stop codons, then translates to amino acid.
//'
//' Uses an inline standard genetic code table (no external codon map needed).
//'
//' @param junction_v_genes Character vector. V gene names for each junction.
//' @param junction_j_genes Character vector. J gene names for each junction.
//' @param junction_nucseqs Character vector. CDR3 nucleotide sequences.
//' @param junction_breakpoints_pre_d List of integer vectors. Pre-D breakpoints
//'   for each junction.
//' @param junction_breakpoints_post_d List of integer vectors. Post-D
//'   breakpoints for each junction.
//' @param chain Character string. \code{"A"} for alpha (pre-D breakpoints only)
//'   or \code{"B"} for beta (both pre-D and post-D).
//' @param num_samples Integer. Number of chimeric TCRs to generate.
//' @param max_attempts Integer. Maximum sampling attempts before returning
//'   partial results.
//' @return A data.frame with columns: \code{v_gene}, \code{j_gene},
//'   \code{cdr3} (amino acid), \code{cdr3_nucseq}.
//' @examples
//' \dontrun{
//'   result <- rcpp_resample_shuffled_tcr_chains(
//'     v_genes, j_genes, nucseqs, bp_pre, bp_post, "B", 1000L, 100000L
//'   )
//' }
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame rcpp_resample_shuffled_tcr_chains(
    const CharacterVector& junction_v_genes,
    const CharacterVector& junction_j_genes,
    const CharacterVector& junction_nucseqs,
    const Rcpp::List& junction_breakpoints_pre_d,
    const Rcpp::List& junction_breakpoints_post_d,
    const std::string& chain,
    int num_samples,
    int max_attempts) {

    int n_junctions = junction_v_genes.size();
    if (n_junctions == 0) {
        Rcpp::stop("rcpp_resample_shuffled_tcr_chains: no junctions provided");
    }
    if (chain != "A" && chain != "B") {
        Rcpp::stop("rcpp_resample_shuffled_tcr_chains: chain must be 'A' or 'B', got '%s'",
                   chain.c_str());
    }

    // Build inline codon translation table
    std::unordered_map<std::string, char> codon_table = init_codon_table();

    // Pre-convert to C++ types
    std::vector<std::string> cpp_v_genes(n_junctions);
    std::vector<std::string> cpp_j_genes(n_junctions);
    std::vector<std::string> cpp_nucseqs(n_junctions);
    std::vector<std::vector<int>> cpp_bp_pre(n_junctions);
    std::vector<std::vector<int>> cpp_bp_post(n_junctions);

    for (int i = 0; i < n_junctions; ++i) {
        cpp_v_genes[i] = Rcpp::as<std::string>(junction_v_genes[i]);
        cpp_j_genes[i] = Rcpp::as<std::string>(junction_j_genes[i]);
        cpp_nucseqs[i] = Rcpp::as<std::string>(junction_nucseqs[i]);

        IntegerVector bp_pre = junction_breakpoints_pre_d[i];
        IntegerVector bp_post = junction_breakpoints_post_d[i];
        cpp_bp_pre[i].assign(bp_pre.begin(), bp_pre.end());
        cpp_bp_post[i].assign(bp_post.begin(), bp_post.end());
    }

    // Output storage
    std::vector<std::string> out_v_gene(num_samples);
    std::vector<std::string> out_j_gene(num_samples);
    std::vector<std::string> out_cdr3(num_samples);
    std::vector<std::string> out_nucseq(num_samples);

    int count = 0;
    int attempts = 0;
    int successes = 0;

    while (count < num_samples && attempts < max_attempts) {
        // Sample two junctions using R's RNG (for set.seed reproducibility)
        int idx1 = static_cast<int>(R::runif(0.0, static_cast<double>(n_junctions)));
        int idx2 = static_cast<int>(R::runif(0.0, static_cast<double>(n_junctions)));
        if (idx1 >= n_junctions) idx1 = n_junctions - 1;
        if (idx2 >= n_junctions) idx2 = n_junctions - 1;

        const std::string& nucseq1 = cpp_nucseqs[idx1];
        const std::string& nucseq2 = cpp_nucseqs[idx2];
        if (nucseq1 == nucseq2) { ++attempts; continue; }

        // Choose breakpoint indices to try
        int inds[2];
        int n_inds;
        if (chain == "A") {
            inds[0] = 0;  // only pre_d
            n_inds = 1;
        } else {
            // Randomize order of pre_d vs post_d
            if (R::runif(0.0, 1.0) < 0.5) {
                inds[0] = 0; inds[1] = 1;
            } else {
                inds[0] = 1; inds[1] = 0;
            }
            n_inds = 2;
        }

        bool success = false;
        for (int ii = 0; ii < n_inds; ++ii) {
            const std::vector<int>& bp_set1 =
                (inds[ii] == 0) ? cpp_bp_pre[idx1] : cpp_bp_post[idx1];
            const std::vector<int>& bp_set2 =
                (inds[ii] == 0) ? cpp_bp_pre[idx2] : cpp_bp_post[idx2];

            // Find shared breakpoints
            // Use a set for the smaller list
            std::vector<int> shared;
            if (bp_set1.size() <= bp_set2.size()) {
                std::unordered_set<int> s1(bp_set1.begin(), bp_set1.end());
                for (int bp : bp_set2) {
                    if (s1.count(bp)) shared.push_back(bp);
                }
            } else {
                std::unordered_set<int> s2(bp_set2.begin(), bp_set2.end());
                for (int bp : bp_set1) {
                    if (s2.count(bp)) shared.push_back(bp);
                }
            }

            if (shared.empty()) continue;

            // Pick random shared breakpoint
            int bp_idx = static_cast<int>(R::runif(0.0, static_cast<double>(shared.size())));
            if (bp_idx >= static_cast<int>(shared.size())) bp_idx = static_cast<int>(shared.size()) - 1;
            int bp = shared[bp_idx];

            // Construct chimeric nucleotide sequence
            std::string nucseq;
            int len1 = static_cast<int>(nucseq1.size());
            int len2 = static_cast<int>(nucseq2.size());

            if (bp >= 0) {
                // nucseq1[:bp] + nucseq2[bp:]
                nucseq = nucseq1.substr(0, bp) + nucseq2.substr(bp);
            } else {
                // nucseq1[:len1+bp] + nucseq2[len2+bp:]
                int cut1 = len1 + bp;
                int cut2 = len2 + bp;
                if (cut1 < 0) cut1 = 0;
                if (cut2 < 0) cut2 = 0;
                nucseq = nucseq1.substr(0, cut1) + nucseq2.substr(cut2);
            }

            // Check length divisible by 3
            if (nucseq.size() % 3 != 0) continue;

            // Translate and check for stop codons inline
            int naa = static_cast<int>(nucseq.size()) / 3;
            std::string cdr3;
            cdr3.reserve(naa);
            bool has_stop = false;

            for (int ci = 0; ci < naa; ++ci) {
                std::string codon = nucseq.substr(ci * 3, 3);
                // Codon table uses uppercase; nucseqs may be lowercase
                for (char& c : codon) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                auto it = codon_table.find(codon);
                char aa;
                if (it != codon_table.end()) {
                    aa = it->second;
                } else {
                    aa = 'X';
                }
                if (aa == '*' || aa == 'X') {
                    has_stop = true;
                    break;
                }
                cdr3.push_back(aa);
            }

            if (!has_stop) {
                success = true;
                out_v_gene[count] = cpp_v_genes[idx1];
                out_j_gene[count] = cpp_j_genes[idx2];
                out_cdr3[count]   = cdr3;
                out_nucseq[count] = nucseq;
                ++count;
                break;
            }
        }

        ++attempts;
        if (success) ++successes;
        if (attempts % 1000 == 0) Rcpp::checkUserInterrupt();
        if (success && count >= num_samples) break;
    }

    if (count < num_samples) {
        Rf_warning(
            "rcpp_resample_shuffled_tcr_chains: reached max_attempts (%d) with only "
            "%d/%d samples (success_rate=%.2f%%). Returning partial results.",
            max_attempts, count, num_samples,
            attempts > 0 ? 100.0 * successes / attempts : 0.0);
    }

    Rprintf("rcpp_resample_shuffled_tcr_chains: success_rate: %.2f\n",
            attempts > 0 ? 100.0 * successes / attempts : 0.0);

    CharacterVector out_v(count), out_j(count), out_c(count), out_n(count);
    for (int i = 0; i < count; ++i) {
        out_v[i] = out_v_gene[i];
        out_j[i] = out_j_gene[i];
        out_c[i] = out_cdr3[i];
        out_n[i] = out_nucseq[i];
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("v_gene")      = out_v,
        Rcpp::Named("j_gene")      = out_j,
        Rcpp::Named("cdr3")        = out_c,
        Rcpp::Named("cdr3_nucseq") = out_n,
        Rcpp::Named("stringsAsFactors") = false
    );
}
