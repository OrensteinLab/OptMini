


#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <shared_mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <omp.h>

// --- Project headers ---
#include "config.h"
#include "functions.h"
#include "tools.h"
#include "density_and_gc.h"
#include "Minimizers.h"
#include "common_structures.h"
#include "uhs.h"
#include "dp.h"
#include "dfs.h"
#include "random_dens.h"



std::pair<uint64_t, uint64_t> exhaust(uint32_t w, uint32_t k) {


	std::unordered_map<uint32_t, uint64_t> gc_upper_bound = {
		{3, 105}, {4, 172}, {5, 282}, {6, 483}, {7, 846}, {8, 1514},
		{9, 2741}, {10, 5014}, {11, 9239}, {12, 17066}, {13, 31741},
		{14, 59410}, {15, 111817}
	};



	std::cout << "gm_upper_bound:\t" << gc_upper_bound[w] << std::endl;

	uint64_t n_kmers = static_cast<uint64_t>(1) << k;
	uint64_t kmer_mask = n_kmers - 1;

	std::vector<uint64_t> gc(n_kmers, static_cast<uint64_t>(1) << w);
	for (uint64_t prefix = 0; prefix < n_kmers; prefix++)
	{
		uint64_t next_kmer = (prefix << 1) & kmer_mask;
		dfs_initial_wrapper(next_kmer, { prefix }, 1, w, k, gc);
		dfs_initial_wrapper(next_kmer | 1, { prefix }, 1, w, k, gc);
	}

	// #key is a mask of a subset, value is gc + correction for 0..0 and 1..1
	std::unordered_map<uint64_t, uint64_t> old_mask_to_gc;

	uint64_t n_kmers_starting_with_0 = static_cast<uint64_t>(1) << (k - 1);

	old_mask_to_gc[1] = gc[0] + 1; // only need correction for 11111
	for (uint64_t prefix = 1; prefix < n_kmers_starting_with_0; prefix++)
	{
		old_mask_to_gc[static_cast<uint64_t>(1) << prefix] = gc[prefix] + 2;
	}






	// now for the main cycle
	for (uint64_t rank = 0; rank <= n_kmers; rank++)
	{
		std::unordered_map<uint64_t, uint64_t> mask_to_gc;

		// loop over the old_mask_to_gc
		for (std::pair<const uint64_t, uint64_t>& entry : old_mask_to_gc)
		{
			uint64_t mask = entry.first;
			uint64_t value = entry.second;
			bool can_be_uhs = true;

			for (uint64_t kmer = 0; kmer < n_kmers; kmer++)
			{
				uint64_t temp_value = value;
				// if kmer isn't in the mask
				if (!((mask >> kmer) & 1))
				{
					// handle 1111 or 00000
					if (kmer == 0 || kmer == (n_kmers - 1))
					{
						temp_value -= 1;
					}
					uint64_t new_mask = mask | (static_cast<uint64_t>(1) << kmer);
					uint64_t next_kmer = (kmer << 1) & kmer_mask;
					uint64_t previous_kmer = kmer >> 1;
					uint64_t prefix_to_add = dfs_prefix_wrapper(next_kmer, mask, 1, w, k) + dfs_prefix_wrapper(next_kmer | 1, mask, 1, w, k);
					uint64_t suffix_to_add = 0;
					temp_value += prefix_to_add;
					// todo: might need ot be <=
					if (temp_value <= gc_upper_bound[w])
					{
						suffix_to_add = dfs_suffix_wrapper(previous_kmer, new_mask, 1, w, k) + dfs_suffix_wrapper(previous_kmer | (static_cast<uint64_t>(1) << (k - 1)), new_mask, 1, w, k);
						temp_value += suffix_to_add;
						if (temp_value <= gc_upper_bound[w])
						{
							// updathe the value in the mask_to_gc unless it has a value then take the minimum
							if (mask_to_gc.find(new_mask) == mask_to_gc.end())
							{
								mask_to_gc[new_mask] = temp_value;
							}
							else
							{
								mask_to_gc[new_mask] = std::min(mask_to_gc[new_mask], temp_value);
							}
						}
					}
					if (prefix_to_add > 0 || suffix_to_add > 0) {
						can_be_uhs = false;
					}

				}
			}
			if (can_be_uhs)
			{
				// check if the mask is UHS
				if (is_uhs_dfs_wrapper(mask, w, k))
				{
					std::cout << "UHS found" << std::endl;
					std::cout << "mask:\t" << std::bitset<64>(mask) << "\tgc:\t" << value << std::endl;
					return std::make_pair(mask, value);
				}
			}

		}

		uint32_t set_size = rank + 2;
		old_mask_to_gc = mask_to_gc;
		uint64_t n_masks = old_mask_to_gc.size();



		uint64_t max_gc = std::max_element(old_mask_to_gc.begin(), old_mask_to_gc.end(),
			[](const auto& a, const auto& b) { return a.second < b.second; })->second;

		uint64_t min_gc = std::min_element(old_mask_to_gc.begin(), old_mask_to_gc.end(),
			[](const auto& a, const auto& b) { return a.second < b.second; })->second;

		std::cout << "mask size\t" << set_size << "\tn_masks:\t" << old_mask_to_gc.size() << "\tmin_gc:\t" << min_gc << "\tmax_gc:\t" << max_gc << std::endl;


		if (n_masks == 0)
		{
			std::cout << "upper bound is optimal" << std::endl;
			return std::make_pair(0, 0);
		}

	}

	return std::make_pair(0, 0);
}



