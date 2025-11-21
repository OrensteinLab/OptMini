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


// -------------------------------------------------------------
// DFS utilities
// -------------------------------------------------------------

// Initial recursive DFS for counting GC contexts
void dfs_initial(uint64_t kmer,
    std::unordered_set<uint64_t> kmers,
    uint32_t levels_left,
    uint64_t kmer_mask,
    std::vector<uint64_t>& gc);

// Wrapper to simplify initial DFS call
void dfs_initial_wrapper(uint64_t kmer,
    std::unordered_set<uint64_t> kmers,
    uint32_t lev,
    uint32_t w,
    uint32_t k,
    std::vector<uint64_t>& gc);

// Recursive prefix DFS
uint64_t dfs_prefix(uint64_t kmer,
    uint64_t mask,
    uint32_t levels_left,
    uint64_t kmer_mask);

// Wrapper for prefix DFS
uint64_t dfs_prefix_wrapper(uint64_t kmer,
    uint64_t mask,
    uint32_t lev,
    uint32_t w,
    uint32_t k);

// Recursive suffix DFS
uint64_t dfs_suffix(uint64_t kmer,
    uint64_t mask,
    uint32_t levels_left,
    uint64_t left_bit);

// Wrapper for suffix DFS
uint64_t dfs_suffix_wrapper(uint64_t kmer,
    uint64_t mask,
    uint32_t lev,
    uint32_t w,
    uint32_t k);

// -------------------------------------------------------------
// Exhaustive UHS exploration
// -------------------------------------------------------------

// Single-threaded baseline (original DFS exhaustive version)
std::pair<uint64_t, uint64_t> dij_exhaust_dfs(uint32_t w,
    uint32_t k,
    uint64_t n_charged);

// Multi-threaded shared-heap DFS exploration
std::pair<uint64_t, uint64_t> dij_exhaust_parallel_shared_heap(uint32_t w,
    uint32_t k,
    uint64_t n_charged,
    int N = 8);

// Single-thread fallback for remaining heap processing
static void process_remaining_heap_single_dfs(
    std::vector<std::pair<uint64_t, uint64_t>>& globalHeap,
    std::array<SeenShard, N_SHARDS>& seen_shards,
    uint32_t w,
    uint32_t k,
    uint64_t n_charged,
    std::atomic<uint64_t>& best_gc_found,
    std::atomic<uint64_t>& best_mask_found);

// -------------------------------------------------------------
// Reconstruction utilities
// -------------------------------------------------------------

// Given a final mask + best GC, reconstruct the order of kmers
std::vector<uint64_t> get_best_order(uint64_t original_mask,
    uint64_t best_gc,
    uint32_t w,
    uint32_t k);