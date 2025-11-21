#pragma once
#include <iostream>
#include <chrono>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <limits>
#include "functions.h" 
#include "tools.h"
#include "density_and_gc.h"
#include <algorithm> 
#include <cmath>
#include <random>     
#include <map>
#include <iomanip> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "config.h"
#include <iostream>
#include <cstdlib>  
#include <ctime>    
#include <cstdlib>  // For std::stoul, std::stod
#include <stdexcept>
#include "Minimizers.h"
#include <unordered_set>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <omp.h>
#include <iostream>
#include <utility> // for std::pair
#include <iostream>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <vector>
#include <queue>
#include <unordered_map>
#include <utility>
#include <cstdint>
#include <atomic>
#include <array>
#include <functional>
#include <condition_variable>
#include <chrono>
#include "uhs.h"
#include "common_structures.h"


#ifndef DP_H
#define DP_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <array>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>
#include <numeric>
#include <algorithm>
#include <cstdint>
#include <utility>
#include <chrono>

#include "common_structures.h"
#include "uhs.h"
#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::uint256_t;


// =====================================================
// =============== 64-bit DP Functions ==================
// =====================================================

// Core dynamic-programming GC computation
std::pair<uint64_t, uint64_t> dij_exhaust_dp(uint32_t w, uint32_t k, uint64_t n_charged);

// Parallel shared-heap search
std::pair<uint64_t, uint64_t>
dij_exhaust_parallel_shared_heap_dp(uint32_t w, uint32_t k, uint64_t n_charged, int N = 8);

// Parallel shared-heap search (restrict useless kmers)
std::pair<uint64_t, uint64_t>
dij_exhaust_parallel_shared_heap_dp_restrict_masks(uint32_t w, uint32_t k, uint64_t n_charged, int N = 8);

// Reconstruct best order from mask/gc
std::vector<uint64_t> get_best_order_dp(uint64_t original_mask, uint64_t best_gc,
    uint32_t w, uint32_t k);

static void process_remaining_heap_single_dp(
	std::vector<std::pair<uint64_t, uint64_t>>& globalHeap,
	std::array<SeenShard, N_SHARDS>& seen_shards,
	uint32_t w, uint32_t k, uint64_t n_charged,
	std::atomic<uint64_t>& best_gc_found,
	std::atomic<uint64_t>& best_mask_found);

static void process_remaining_heap_single_dp_restrict_masks(
	std::vector<std::pair<uint64_t, uint64_t>>& globalHeap,
	std::array<SeenShard, N_SHARDS>& seen_shards,
	uint32_t w, uint32_t k, uint64_t n_charged,
	std::atomic<uint64_t>& best_gc_found,
	std::atomic<uint64_t>& best_mask_found);

// =====================================================
// =============== 256-bit DP Functions =================
// =====================================================

uint256_t prefix_charged_dp_256(uint64_t mask, uint64_t kmer_to_add, uint32_t w, uint32_t k);

uint256_t suffix_charged_dp_256_v2(uint64_t old_mask, uint64_t kmer_to_add,
	uint32_t w, uint32_t k);

std::pair<uint64_t, boost::multiprecision::uint256_t> dij_exhaust_dp_256_time_limit(uint32_t w, uint32_t k, boost::multiprecision::uint256_t n_charged, double time_limit_minutes, uint32_t clean_every_n);

// Sequential exhaustive DP for 256-bit masks
std::pair<uint64_t, uint256_t>
dij_exhaust_dp_256(uint32_t w, uint32_t k, uint256_t n_charged);



//// Parallel shared-heap (256-bit)
//std::pair<uint64_t, uint256_t>
//dij_exhaust_parallel_shared_heap_dp_256(uint32_t w, uint32_t k, uint256_t n_charged, int N = 8);
//
//// Parallel shared-heap with restricted masks (256-bit)
//std::pair<uint64_t, uint256_t>
//dij_exhaust_parallel_shared_heap_dp_restrict_masks_256(uint32_t w, uint32_t k,
//    uint256_t n_charged, int N = 8);
	

// Reconstruct best order for 256-bit version
std::vector<uint64_t> get_best_order_dp_256(uint64_t original_mask,
	boost::multiprecision::uint256_t best_gc,
	uint32_t w, uint32_t k);

std::pair<std::vector<uint64_t>, uint256_t> get_best_order_and_gc_dp_256(uint64_t original_mask, boost::multiprecision::uint256_t best_gc, uint32_t w, uint32_t k, bool verbose);



//static void process_remaining_heap_single_dp_256(
//	std::vector<std::pair<uint256_t, uint64_t>>& heap,
//	std::array<SeenShard64x256, N_SHARDS>& seen_shards,
//	uint32_t w, uint32_t k, uint256_t n_charged,
//	std::atomic<uint256_t>& best_gc,
//	std::atomic<uint64_t>& best_mask);
//
//static void process_remaining_heap_single_dp_restrict_masks_256(
//	std::vector<std::pair<uint256_t, uint64_t>>& heap,
//	std::array<SeenShard64x256, N_SHARDS>& seen_shards,
//	uint32_t w, uint32_t k, uint256_t n_charged,
//	std::atomic<uint256_t>& best_gc,
//	std::atomic<uint64_t>& best_mask);
	

#endif // DP_H

