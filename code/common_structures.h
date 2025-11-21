#pragma once
// Custom comparator for min-heap (based on gc value)

#include <iostream>
#include <chrono>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <limits>
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
#include <boost/multiprecision/cpp_dec_float.hpp>
//#include "functions.h" 




struct MinHeapComparator {
	bool operator()(const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) {
		return a.first > b.first; // Compare gc values for min-heap
	}
};

// For 256-bit unsigned integers
using boost::multiprecision::uint256_t;

struct MinHeapComparator_256 {
	bool operator()(const std::pair<uint256_t, uint256_t>& a,
		const std::pair<uint256_t, uint256_t>& b) const noexcept {
		return a.first > b.first; // Compare gc values for min-heap
	}
};




// ===================== CONFIG =====================
constexpr size_t N_SHARDS = 100;
constexpr size_t WORKLOAD_PER_WORKER = 10000;
constexpr size_t MAX_LOCAL_PUSH_BATCH = 1000;
constexpr int    CLEANUP_PERIOD_MS = 2000;
constexpr bool   VERBOSE_STATUS = true;


// ===================== STRIPED SEEN_MASKS =====================
struct SeenShard {
	std::unordered_map<uint64_t, uint64_t> map;
	std::mutex lock;
};

static inline size_t shard_id_for(uint64_t key) {
	return std::hash<uint64_t>{}(key) % N_SHARDS;
}

static inline bool seen_lookup(std::array<SeenShard, N_SHARDS>& shards,
	uint64_t key, uint64_t& out_val)
{
	const size_t sid = shard_id_for(key);
	auto& sh = shards[sid];
	std::lock_guard<std::mutex> g(sh.lock);
	auto it = sh.map.find(key);
	if (it == sh.map.end()) return false;
	out_val = it->second;
	return true;
}

static inline bool seen_update_if_better(std::array<SeenShard, N_SHARDS>& shards,
	uint64_t key, uint64_t val)
{
	const size_t sid = shard_id_for(key);
	auto& sh = shards[sid];
	std::lock_guard<std::mutex> g(sh.lock);
	auto it = sh.map.find(key);
	if (it == sh.map.end() || val < it->second) {
		sh.map[key] = val;
		return true;
	}
	return false;
}

static inline void seen_set(std::array<SeenShard, N_SHARDS>& shards,
	uint64_t key, uint64_t val)
{
	const size_t sid = shard_id_for(key);
	auto& sh = shards[sid];
	std::lock_guard<std::mutex> g(sh.lock);
	sh.map[key] = val;
}


// ===================== REBUILD HEAP =====================
static inline void rebuild_minheap_dedup(std::vector<std::pair<uint64_t, uint64_t>>& heap,
	std::array<SeenShard, N_SHARDS>& seen_shards)
{
	std::unordered_map<uint64_t, uint64_t> best;
	best.reserve(heap.size() * 1.3);

	for (auto& [val, m] : heap) {
		auto it = best.find(m);
		if (it == best.end() || val < it->second) best[m] = val;
	}

	heap.clear();
	heap.reserve(best.size());
	for (auto& [m, v] : best) {
		heap.emplace_back(v, m);
		seen_update_if_better(seen_shards, m, v);
	}
	std::make_heap(heap.begin(), heap.end(), MinHeapComparator{});
}

// ===================== UHS RESULT =====================
struct UhsCandidate {
	bool     found = false;
	uint64_t mask = 0;
	uint64_t gc = 0;
};






// Custom hash for uint256_t (since std::hash doesn’t support it)
struct UInt256Hash {
	size_t operator()(const uint256_t& x) const noexcept {
		// Use lower 64 bits for hash
		return std::hash<uint64_t>{}(x.convert_to<uint64_t>());
	}
};

// Equality operator (required by unordered_map)
struct UInt256Eq {
	bool operator()(const uint256_t& a, const uint256_t& b) const noexcept {
		return a == b;
	}
};

struct SeenShard_256 {
	std::unordered_map<uint256_t, uint256_t, UInt256Hash, UInt256Eq> map;
	std::mutex lock;
};

static inline size_t shard_id_for_256(const uint256_t& key) {
	return UInt256Hash{}(key) % N_SHARDS;
}

static inline bool seen_lookup_256(std::array<SeenShard_256, N_SHARDS>& shards,
	const uint256_t& key, uint256_t& out_val)
{
	const size_t sid = shard_id_for_256(key);
	auto& sh = shards[sid];
	std::lock_guard<std::mutex> g(sh.lock);
	auto it = sh.map.find(key);
	if (it == sh.map.end()) return false;
	out_val = it->second;
	return true;
}

static inline bool seen_update_if_better_256(std::array<SeenShard_256, N_SHARDS>& shards,
	const uint256_t& key, const uint256_t& val)
{
	const size_t sid = shard_id_for_256(key);
	auto& sh = shards[sid];
	std::lock_guard<std::mutex> g(sh.lock);
	auto it = sh.map.find(key);
	if (it == sh.map.end() || val < it->second) {
		sh.map[key] = val;
		return true;
	}
	return false;
}

static inline void seen_set_256(std::array<SeenShard_256, N_SHARDS>& shards,
	const uint256_t& key, const uint256_t& val)
{
	const size_t sid = shard_id_for_256(key);
	auto& sh = shards[sid];
	std::lock_guard<std::mutex> g(sh.lock);
	sh.map[key] = val;
}

// ===================== REBUILD HEAP (256-bit) =====================


static inline void rebuild_minheap_dedup_256(
	std::vector<std::pair<uint256_t, uint256_t>>& heap,
	std::array<SeenShard_256, N_SHARDS>& seen_shards)
{
	std::unordered_map<uint256_t, uint256_t, UInt256Hash, UInt256Eq> best;
	best.reserve(heap.size() * 1.3);

	for (auto& [val, m] : heap) {
		auto it = best.find(m);
		if (it == best.end() || val < it->second)
			best[m] = val;
	}

	heap.clear();
	heap.reserve(best.size());
	for (auto& [m, v] : best) {
		heap.emplace_back(v, m);
		seen_update_if_better_256(seen_shards, m, v);
	}
	std::make_heap(heap.begin(), heap.end(), MinHeapComparator_256{});
}

// ===================== UHS RESULT (256-bit) =====================
struct UhsCandidate_256 {
	bool found = false;
	uint256_t mask = 0;
	uint256_t gc = 0;
};




// --- comparator for (gc=uint256_t, mask=uint64_t) min-heap ---
struct MinHeapComparator256_64 {
	bool operator()(const std::pair<boost::multiprecision::uint256_t, uint64_t>& a,
		const std::pair<boost::multiprecision::uint256_t, uint64_t>& b) const noexcept {
		return a.first > b.first; // min-heap by gc
	}
};

struct MinHeapComparator256_64_64 {
	bool operator()(
		const std::pair<boost::multiprecision::uint256_t,
		std::pair<uint64_t, uint64_t>>&a,
		const std::pair<boost::multiprecision::uint256_t,
		std::pair<uint64_t, uint64_t>>&b) const noexcept
	{
		return a.first > b.first;   // min-heap by gc
	}
};



// --- sharded seen map: key=uint64_t mask, val=uint256_t gc ---
struct SeenShard64x256 {
	std::unordered_map<uint64_t, boost::multiprecision::uint256_t> map;
	std::mutex lock;
};

static inline size_t shard_id_for_64(uint64_t key) {
	return std::hash<uint64_t>{}(key) % N_SHARDS;
}

static inline bool seen_lookup_64x256(std::array<SeenShard64x256, N_SHARDS>& shards,
	uint64_t key, boost::multiprecision::uint256_t& out_val)
{
	auto& sh = shards[shard_id_for_64(key)];
	std::lock_guard<std::mutex> g(sh.lock);
	auto it = sh.map.find(key);
	if (it == sh.map.end()) return false;
	out_val = it->second;
	return true;
}

static inline bool seen_update_if_better_64x256(std::array<SeenShard64x256, N_SHARDS>& shards,
	uint64_t key, const boost::multiprecision::uint256_t& val)
{
	auto& sh = shards[shard_id_for_64(key)];
	std::lock_guard<std::mutex> g(sh.lock);
	auto it = sh.map.find(key);
	if (it == sh.map.end() || val < it->second) { sh.map[key] = val; return true; }
	return false;
}

static inline void seen_set_64x256(std::array<SeenShard64x256, N_SHARDS>& shards,
	uint64_t key, const boost::multiprecision::uint256_t& val)
{
	auto& sh = shards[shard_id_for_64(key)];
	std::lock_guard<std::mutex> g(sh.lock);
	sh.map[key] = val;
}

static inline void rebuild_minheap_dedup_256_64(
	std::vector<std::pair<boost::multiprecision::uint256_t, uint64_t>>& heap,
	std::array<SeenShard64x256, N_SHARDS>& seen_shards)
{
	using boost::multiprecision::uint256_t;
	std::unordered_map<uint64_t, uint256_t> best;
	best.reserve(heap.size() * 2);

	for (auto& p : heap) {
		auto& [gc, m] = p;
		auto it = best.find(m);
		if (it == best.end() || gc < it->second) best[m] = gc;
	}

	heap.clear();
	heap.reserve(best.size());
	for (auto& kv : best) {
		heap.emplace_back(kv.second, kv.first);
		seen_update_if_better_64x256(seen_shards, kv.first, kv.second);
	}
	std::make_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
}
