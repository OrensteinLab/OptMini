#pragma once


#include <boost/multiprecision/cpp_int.hpp>
#include <vector>
#include <string>
#include <tuple>
#include "window_types_manager.h"

uint64_t prob_gc(uint32_t w, uint32_t k, const std::vector<uint64_t>& order);
uint64_t gc_dfs(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev);
uint64_t gc_dfs_prev_could_be_min(uint32_t w, uint64_t n_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev);
uint64_t gc_dfs_pref_is_not_min(uint32_t w, uint64_t n_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev);

using namespace boost::multiprecision;
uint64_t prob_gc(uint32_t w, uint32_t k, const std::vector<uint64_t>& order);


double density_expected_binary(uint32_t w, uint32_t k, const std::vector<uint64_t>& order);
