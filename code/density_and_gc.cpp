#include <vector>
#include <cstdint>
#include <iostream>
#include <chrono>
#include <bitset>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>
#include "tools.h"
#include "density_and_gc.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <sstream>  // Include this for istringstream
#include <iostream>
#include <fstream>
#include "functions.h"
#include "config.h"
#include "density_and_gc.h"

using namespace boost::multiprecision;



uint64_t prob_gc(uint32_t w, uint32_t k, const std::vector<uint64_t>& order) {
    uint64_t total = 0;
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;

    for (uint64_t pref = 0; pref < n_kmers; ++pref) {
        uint64_t nextkmer = (pref << 1) & n_kmers_mask;
        total += gc_dfs(w, k, order, nextkmer, order[pref], 1) +
            gc_dfs(w, k, order, nextkmer + 1, order[pref], 1);
    }

    return total;
}

uint64_t gc_dfs(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev) {
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;
    return gc_dfs_prev_could_be_min(w, n_kmers_mask, order, kmer, minorder, lev);
}

uint64_t gc_dfs_prev_could_be_min(uint32_t w, uint64_t n_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev) {
    if (lev == w) {
        return 1;
    }
    else {
        uint64_t nextkm = (kmer << 1) & n_kmers_mask;
        if (order[kmer] >= minorder) {
            return gc_dfs_prev_could_be_min(w, n_kmers_mask, order, nextkm, minorder, lev + 1) +
                gc_dfs_prev_could_be_min(w, n_kmers_mask, order, nextkm + 1, minorder, lev + 1);
        }
        else {
            return gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm, order[kmer], lev + 1) +
                gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm + 1, order[kmer], lev + 1);
        }
    }
}

uint64_t gc_dfs_pref_is_not_min(uint32_t w, uint64_t n_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev) {
    if (lev == w) {
        return order[kmer] < minorder;
    }
    else {
        uint64_t nextkm = (kmer << 1) & n_kmers_mask;
        if (order[kmer] >= minorder) {
            return gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm, minorder, lev + 1) +
                gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm + 1, minorder, lev + 1);
        }
        else {
            return gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm, order[kmer], lev + 1) +
                gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm + 1, order[kmer], lev + 1);
        }
    }
}

double density_expected_binary(uint32_t w, uint32_t k, const std::vector<uint64_t>& order) {
	uint64_t total = prob_gc(w, k, order);
	double density = calc_density(total, k, w);
	return density;
}
