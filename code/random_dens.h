#pragma once
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


#include <iostream>
#include <thread>
#include <mutex>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cstdint>
#include <atomic>
#include <array>
#include <functional>
#include <condition_variable>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>


using boost::multiprecision::uint256_t;
using boost::multiprecision::cpp_dec_float_100;


// ---------------------------------------------------------
// Function declarations
// ---------------------------------------------------------

// === 64-bit versions ===

// Compute density curve for given (k, w)
std::vector<double> randens_w(uint32_t k, uint32_t w);

// Return the final random density value (for w-2)
double get_random_density(uint32_t w, uint32_t k);


// === 256-bit versions ===

// Compute density curve using 256-bit integer arithmetic
std::vector<double> randens_w_256(uint32_t k, uint32_t w);

// Return the final random density (256-bit internal version)
double get_random_density_256(uint32_t w, uint32_t k);