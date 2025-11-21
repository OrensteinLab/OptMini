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
#include <boost/multiprecision/cpp_int.hpp>
#include "dfs.h"
#include "dp.h"

using boost::multiprecision::uint256_t;



// --- 64-bit versions ---

std::pair<std::list<uint64_t>, uint64_t> get_opt_uhs_dfs(uint32_t w, uint32_t k);
std::pair<std::list<uint64_t>, uint64_t> get_opt_uhs_dp(uint32_t w, uint32_t k);
std::pair<std::list<uint64_t>, uint64_t> get_opt_uhs_dp_256(uint32_t w, uint32_t k);
std::pair<std::vector<uint64_t>, uint64_t> get_best_order_from_uhs(uint64_t original_mask, uint64_t uhs_size, uint32_t w, uint32_t k);

bool is_uhs_dfs(uint64_t kmer, uint64_t mask, uint32_t levels_left, uint64_t kmer_mask);
bool is_uhs_dfs_wrapper(uint64_t mask, uint32_t w, uint32_t k);

bool is_uhs_dp(uint64_t mask, uint32_t w, uint32_t k);
bool is_uhs_dp_wrapper(uint64_t mask, uint32_t w, uint32_t k);

uint64_t how_many_contexts_uncovered_by_mask_dp(uint32_t w, uint32_t k, uint64_t mask);
bool is_kmer_useless_dp(uint32_t w, uint32_t k, uint64_t mask, uint64_t kmer_to_add);


// --- 256-bit versions ---


bool is_uhs_dp_256(uint64_t mask, uint32_t w, uint32_t k);
bool is_uhs_dp_wrapper_256(uint64_t mask, uint32_t w, uint32_t k);

uint32_t get_uhs_size_dp_256(std::vector<uint64_t> order, uint32_t w, uint32_t k);

uint256_t how_many_contexts_uncovered_by_mask_dp_256(uint32_t w, uint32_t k, uint64_t mask);
bool is_kmer_useless_dp_256(uint32_t w, uint32_t k, uint64_t mask, uint64_t kmer_to_add);