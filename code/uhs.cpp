
#include "uhs.h"

using boost::multiprecision::uint256_t;

bool is_uhs_dfs(uint64_t kmer, uint64_t mask, uint32_t levels_left, uint64_t kmer_mask) {
	// if kmer is in the mask
	if (((mask >> kmer) & 1)) {
		return true;
	}

	if (levels_left == 0) {
		return false;
	}

	uint64_t next_kmer = (kmer << 1) & kmer_mask;
	return is_uhs_dfs(next_kmer, mask, levels_left - 1, kmer_mask) && is_uhs_dfs(next_kmer | 1, mask, levels_left - 1, kmer_mask);
}


bool is_uhs_dfs_wrapper(uint64_t mask, uint32_t w, uint32_t k) {
	uint64_t kmer_mask = (static_cast<uint64_t>(1) << k) - 1;
	uint64_t n_kmers = static_cast<uint64_t>(1) << k;
	uint32_t levels_left = w - 1;

	for (uint64_t kmer = 0; kmer < n_kmers; kmer++)
	{
		if (!is_uhs_dfs(kmer, mask, levels_left, kmer_mask)) {
			return false;
		}

	}
	return true;
}



bool is_uhs_dp(uint64_t mask, uint32_t w, uint32_t k) {
	uint64_t n_kmers = 1ULL << k;
	uint64_t kmer_mask = n_kmers - 1ULL;

	// we count how many contexts can be constructed using only kmers not in the mask
	std::vector<uint64_t> old_arr(n_kmers, 0ULL);
	std::vector<uint64_t> new_arr(n_kmers, 0ULL);

	// set 1 to all kmers not in the mask
	for (uint64_t kmer = 0; kmer < n_kmers; kmer++) {
		if (((mask >> kmer) & 1) == 0) {
			old_arr[kmer] = 1ULL;
		}
	}

	// run for w-1 times - as we only need to hit every w+k-1 long windows and not a context
	for (uint32_t i = 1; i < w; i++) {
		for (uint64_t kmer = 0; kmer < n_kmers; kmer++) {
			// skip kmers in the mask
			if (((mask >> kmer) & 1)) {
				new_arr[kmer] = 0ULL;
			}
			else {
				uint64_t next_kmer0 = (kmer << 1) & kmer_mask;
				uint64_t next_kmer1 = next_kmer0 | 1ULL;
				new_arr[kmer] = old_arr[next_kmer0] + old_arr[next_kmer1];
			}

		}
		std::swap(old_arr, new_arr);
		// no need to rest to 0s since we overwrite all values
	}

	// check the sum of the array
	uint64_t total = 0ULL;
	for (uint64_t kmer = 0; kmer < n_kmers; kmer++) {
		total += old_arr[kmer];
	}

	// if the sum is positive, return false
	if (total > 0ULL) {
		return false;
	}
	else {
		return true;
	}
}


bool is_uhs_dp_wrapper(uint64_t mask, uint32_t w, uint32_t k) {
	return is_uhs_dp(mask, w, k);
}

using boost::multiprecision::uint256_t;

bool is_uhs_dp_256(uint64_t mask, uint32_t w, uint32_t k) {
	uint64_t n_kmers = 1ULL << k;         // still 64-bit, since k ≤ 64
	uint64_t kmer_mask = n_kmers - 1ULL;

	// we count how many contexts can be constructed using only kmers not in the mask
	std::vector<uint256_t> old_arr(n_kmers, 0);
	std::vector<uint256_t> new_arr(n_kmers, 0);

	// set 1 to all kmers not in the mask
	for (uint64_t kmer = 0; kmer < n_kmers; kmer++) {
		if (((mask >> kmer) & 1) == 0) {
			old_arr[kmer] = 1;
		}
	}

	// run for w-1 times – as we only need to hit every w+k-1-long window
	for (uint32_t i = 1; i < w; i++) {
		for (uint64_t kmer = 0; kmer < n_kmers; kmer++) {
			if (((mask >> kmer) & 1)) {
				new_arr[kmer] = 0;
			}
			else {
				uint64_t next_kmer0 = (kmer << 1) & kmer_mask;
				uint64_t next_kmer1 = next_kmer0 | 1ULL;
				new_arr[kmer] = old_arr[next_kmer0] + old_arr[next_kmer1];
			}
		}
		std::swap(old_arr, new_arr);
	}

	// check the sum of the array
	uint256_t total = 0;
	for (uint64_t kmer = 0; kmer < n_kmers; kmer++) {
		total += old_arr[kmer];
	}



	return total == 0;  // true if no uncovered contexts remain
}

bool is_uhs_dp_wrapper_256(uint64_t mask, uint32_t w, uint32_t k) {
	//std::cout << "trying to check UHS" << std::endl;
	return is_uhs_dp_256(mask, w, k);
}



uint32_t get_uhs_size_dp_256(std::vector<uint64_t> order, uint32_t w, uint32_t k)
{
	uint64_t mask = 0;
	uint64_t uhs_size = 0;
	for (uint64_t kmer : order) {
		// get the added gc
		const uint256_t pre = prefix_charged_dp_256(mask, kmer, w, k);
		const uint256_t suf = suffix_charged_dp_256_v2(mask, kmer, w, k);
		if (pre + suf > 0) {
			uhs_size++;
		}
		mask |= (1ULL << kmer);
	}

	return uhs_size;
}










uint64_t how_many_contexts_uncovered_by_mask_dp(uint32_t w, uint32_t k, uint64_t mask)
{
	uint64_t n_kmers = static_cast<uint64_t>(1) << k;
	uint64_t kmer_mask = n_kmers - 1;
	uint64_t total = 0;

	// array means: how many contexts ending in this kmer are uncovered
	std::vector<uint64_t> old_arr(n_kmers, static_cast<uint64_t>(0));
	std::vector<uint64_t> new_arr(n_kmers, static_cast<uint64_t>(0));

	// populate old_arr with 1 if kmer is not in mask
	for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
		if (!((mask >> kmer) & 1ULL)) {
			old_arr[kmer] = 1ULL;
		}
	}

	// run for w times (since we already did w=0 and need w+1 kmers)
	for (uint64_t i = 0; i < w; ++i) {
		for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
			// if kmer is not in mask
			if (!((mask >> kmer) & 1ULL)) {
				uint64_t prev_kmer = kmer >> 1;
				// sum the previous ones
				new_arr[kmer] = old_arr[prev_kmer] + old_arr[prev_kmer | (1ULL << (k - 1))];
			}
		}
		std::swap(old_arr, new_arr);
		std::fill(new_arr.begin(), new_arr.end(), 0ULL);
	}

	// total should be the sum of all old_arr entries
	total = std::accumulate(old_arr.begin(), old_arr.end(), static_cast<uint64_t>(0));

	return total;
}


// if here we can assume it doesnt charge anything so we need to check only if it exists in the middle of a window where no kmers from mask exist (in at least 1 window)
// assume w >= 2
bool is_kmer_useless_dp(uint32_t w, uint32_t k, uint64_t mask, uint64_t kmer_to_add)
{
	uint64_t n_kmers = 1ULL << k;
	uint64_t kmer_mask = n_kmers - 1ULL;
	uint64_t new_mask = mask | (1ULL << kmer_to_add);

	uint64_t uncovered_1 = how_many_contexts_uncovered_by_mask_dp(w, k, mask);
	uint64_t uncovered_2 = how_many_contexts_uncovered_by_mask_dp(w, k, new_mask);

	if (uncovered_1 == uncovered_2)
	{
		return true;
	}
	else
	{
		return false;
	}

}










// =====================================================
// 256-bit version of how_many_contexts_uncovered_by_mask_dp
// =====================================================
uint256_t how_many_contexts_uncovered_by_mask_dp_256(uint32_t w, uint32_t k, uint64_t mask)
{
	uint64_t n_kmers = 1ULL << k;
	uint64_t kmer_mask = n_kmers - 1ULL;
	uint256_t total = 0;

	std::vector<uint256_t> old_arr(n_kmers, 0);
	std::vector<uint256_t> new_arr(n_kmers, 0);

	// Initialize uncovered kmers
	for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
		if (!((mask >> kmer) & 1)) {
			old_arr[kmer] = 1;
		}
	}

	// Dynamic programming for w steps
	for (uint32_t i = 0; i < w; ++i) {
		for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
			if (!((mask >> kmer) & 1)) {
				uint64_t prev_kmer = kmer >> 1;
				new_arr[kmer] = old_arr[prev_kmer] + old_arr[prev_kmer | (1ULL << (k - 1))];
			}
		}
		std::swap(old_arr, new_arr);
		std::fill(new_arr.begin(), new_arr.end(), 0);
	}

	// Sum uncovered contexts
	for (uint64_t kmer = 0; kmer < n_kmers; ++kmer)
		total += old_arr[kmer];

	return total;
}

// =====================================================
// 256-bit version of is_kmer_useless_dp
// =====================================================
bool is_kmer_useless_dp_256(uint32_t w, uint32_t k, uint64_t mask, uint64_t kmer_to_add)
{
	uint64_t new_mask = mask | (uint64_t(1) << kmer_to_add);

	// needs to be w-1 and not w to go over (w+k-1) windows and not (w+k)

	uint256_t uncovered_1 = how_many_contexts_uncovered_by_mask_dp_256(w-1, k, mask);
	uint256_t uncovered_2 = how_many_contexts_uncovered_by_mask_dp_256(w-1, k, new_mask);

	return (uncovered_1 == uncovered_2);
}









std::pair<std::list<uint64_t>, uint64_t> get_opt_uhs_dfs(uint32_t w, uint32_t k) {
	// make sure k<=6 or otherwise error
	if (k > 6) {
		std::cout << "k must be less than or equal to 6" << std::endl;
		return { {}, 0 };
	}
	// create an empty list
	std::list<uint64_t> optimal_uhs = std::list<uint64_t>();


	// create all kmers of size k
	uint64_t n_kmers = static_cast<uint64_t>(1) << k;
	std::list<uint64_t> masks = std::list<uint64_t>();


	// craete of all masks containing 1 kmer
	for (uint64_t i = 0; i < n_kmers; i++) {
		masks.push_back(static_cast<uint64_t>(1) << i);
	}

	uint64_t opt_uhs_size = 0;

	// main loop- over mask size
	for (uint32_t mask_size = 1; mask_size <= n_kmers; mask_size++) {
		std::cout << "Checking masks of size:\t" << mask_size << "\tNumber of masks to check:\t" << masks.size() << std::endl;
		std::list<uint64_t> new_masks = std::list<uint64_t>();
		// loop over all masks of size mask_size-1
		for (uint64_t mask : masks) {
			if (is_uhs_dfs_wrapper(mask, w, k)) {
				optimal_uhs.push_back(mask);
			}
			else {
				// loop over all kmers
				for (uint64_t kmer = 0; kmer < n_kmers; kmer++) {
					// if kmer is not in the mask
					if (!((mask >> kmer) & 1)) {
						uint64_t new_mask = mask | (static_cast<uint64_t>(1) << kmer);
						new_masks.push_back(new_mask);
					}
				}

			}
		}
		masks = new_masks;

		// if found at least one uhs, break
		if (optimal_uhs.size() > 0) {
			opt_uhs_size = mask_size;
			break;
		}
	}

	// print how many uhs found
	std::cout << "Found " << optimal_uhs.size() << " optimal uhs for w:\t" << w << "\tk:\t" << k << std::endl;

	return { optimal_uhs, opt_uhs_size };


}




std::pair<std::list<uint64_t>, uint64_t> get_opt_uhs_dp(uint32_t w, uint32_t k) {
	// make sure k<=6 or otherwise error
	if (k > 6) {
		std::cout << "k must be less than or equal to 6" << std::endl;
		return { {}, 0 };
	}
	// create an empty list
	std::list<uint64_t> optimal_uhs = std::list<uint64_t>();


	// create all kmers of size k
	uint64_t n_kmers = static_cast<uint64_t>(1) << k;
	std::list<uint64_t> masks = std::list<uint64_t>();


	// craete of all masks containing 1 kmer
	for (uint64_t i = 0; i < n_kmers; i++) {
		masks.push_back(static_cast<uint64_t>(1) << i);
	}

	uint64_t opt_uhs_size = 0;

	// main loop- over mask size
	for (uint32_t mask_size = 1; mask_size <= n_kmers; mask_size++) {
		std::cout << "Checking masks of size:\t" << mask_size << "\tNumber of masks to check:\t" << masks.size() << std::endl;
		std::list<uint64_t> new_masks = std::list<uint64_t>();
		for (uint64_t mask : masks) {
			if (is_uhs_dp_wrapper(mask, w, k)) {
				optimal_uhs.push_back(mask);
			}
			else {
				// loop over all kmers
				for (uint64_t kmer = 0; kmer < n_kmers; kmer++) {
					// if kmer is not in the mask
					if (!((mask >> kmer) & 1)) {
						uint64_t new_mask = mask | (static_cast<uint64_t>(1) << kmer);
						new_masks.push_back(new_mask);
					}
				}

			}
		}
		masks = new_masks;

		// if found at least one uhs, break
		if (optimal_uhs.size() > 0) {
			opt_uhs_size = mask_size;
			break;
		}
	}

	// print how many uhs found
	std::cout << "Found " << optimal_uhs.size() << " optimal uhs for w:\t" << w << "\tk:\t" << k << std::endl;

	return { optimal_uhs, opt_uhs_size };


}

std::pair<std::list<uint64_t>, uint64_t> get_opt_uhs_dp_256(uint32_t w, uint32_t k) {
	// make sure k<=6 or otherwise error
	if (k > 6) {
		std::cout << "k must be less than or equal to 6" << std::endl;
		return { {}, 0 };
	}
	// create an empty list
	std::list<uint64_t> optimal_uhs = std::list<uint64_t>();


	// create all kmers of size k
	uint64_t n_kmers = static_cast<uint64_t>(1) << k;
	std::list<uint64_t> masks = std::list<uint64_t>();


	// craete of all masks containing 1 kmer
	for (uint64_t i = 0; i < n_kmers; i++) {
		masks.push_back(static_cast<uint64_t>(1) << i);
	}

	uint64_t opt_uhs_size = 0;

	// main loop- over mask size
	for (uint32_t mask_size = 1; mask_size <= n_kmers; mask_size++) {
		std::cout << "Checking masks of size:\t" << mask_size << "\tNumber of masks to check:\t" << masks.size() << std::endl;
		std::list<uint64_t> new_masks = std::list<uint64_t>();
		for (uint64_t mask : masks) {
			if (is_uhs_dp_wrapper_256(mask, w, k)) {
				optimal_uhs.push_back(mask);
			}
			else {
				// loop over all kmers
				for (uint64_t kmer = 0; kmer < n_kmers; kmer++) {
					// if kmer is not in the mask
					if (!((mask >> kmer) & 1)) {
						if (is_kmer_useless_dp_256(w, k, mask, kmer)) {
							continue;
						}
						uint64_t new_mask = mask | (static_cast<uint64_t>(1) << kmer);
						new_masks.push_back(new_mask);
					}
				}

			}
		}
		masks = new_masks;

		// if found at least one uhs, break
		if (optimal_uhs.size() > 0) {
			opt_uhs_size = mask_size;
			break;
		}
	}

	// print how many uhs found
	std::cout << "Found " << optimal_uhs.size() << " optimal uhs for w:\t" << w << "\tk:\t" << k << std::endl;

	return { optimal_uhs, opt_uhs_size };


}

std::pair<std::vector<uint64_t>, uint64_t> get_best_order_from_uhs(uint64_t original_mask, uint64_t uhs_size, uint32_t w, uint32_t k) {

	// show the original mask and uhs size
	//std::cout << "Original mask:\t" << std::bitset<64>(original_mask) << "\tUHS size:\t" << uhs_size << std::endl;


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
	std::unordered_map<uint64_t, std::vector<uint64_t>> mask_to_best_order;

	//uint64_t n_kmers_starting_with_0 = static_cast<uint64_t>(1) << (k - 1);

	if ((original_mask >> 0) & 1) {
		old_mask_to_gc[1] = gc[0] + 1; // only need correction for 11111
		mask_to_best_order[1] = { 0 };
	}

	for (uint64_t prefix = 1; prefix < n_kmers; prefix++) // used to be up to n_kmers_starting_with_0
	{
		if ((original_mask >> prefix) & 1) {
			old_mask_to_gc[static_cast<uint64_t>(1) << prefix] = gc[prefix] + 2;
			mask_to_best_order[static_cast<uint64_t>(1) << prefix] = { prefix };
		}

	}

	std::unordered_map<uint64_t, uint64_t> mask_to_gc;

	// now for the main cycle
	for (uint64_t rank = 0; rank < uhs_size - 1; rank++)
	{

		// loop over the old_mask_to_gc
		for (std::pair<const uint64_t, uint64_t>& entry : old_mask_to_gc)
		{
			uint64_t mask = entry.first;
			uint64_t value = entry.second;


			for (uint64_t kmer = 0; kmer < n_kmers; kmer++)
			{
				if (!((original_mask >> kmer) & 1)) {
					continue;
				}
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
					//std::cout << "old_mask:\t" << std::bitset<64>(mask) << "\tkmer_to_add:\t" << kmer << "\tnew_mask:\t" << std::bitset<64>(new_mask) << std::endl;
					uint64_t next_kmer = (kmer << 1) & kmer_mask;
					uint64_t previous_kmer = kmer >> 1;
					temp_value += dfs_prefix_wrapper(next_kmer, mask, 1, w, k) + dfs_prefix_wrapper(next_kmer | 1, mask, 1, w, k);
					temp_value += dfs_suffix_wrapper(previous_kmer, new_mask, 1, w, k) + dfs_suffix_wrapper(previous_kmer | (static_cast<uint64_t>(1) << (k - 1)), new_mask, 1, w, k);
					// updathe the value in the mask_to_gc unless it has a value then take the minimum
					if (mask_to_gc.find(new_mask) == mask_to_gc.end())
					{
						mask_to_gc[new_mask] = temp_value;
						mask_to_best_order[new_mask] = mask_to_best_order[mask];
						mask_to_best_order[new_mask].push_back(kmer);
					}
					else
					{
						if (temp_value < mask_to_gc[new_mask]) {
							mask_to_gc[new_mask] = temp_value;
							mask_to_best_order[new_mask] = mask_to_best_order[mask];
							mask_to_best_order[new_mask].push_back(kmer);
						}
					}
				}

			}

		}

		uint32_t set_size = rank + 2;
		old_mask_to_gc = mask_to_gc;
		uint64_t n_masks = old_mask_to_gc.size();






	}

	// get the best order from the original mask
	std::vector<uint64_t> best_order = mask_to_best_order[original_mask];
	uint64_t best_gc = mask_to_gc[original_mask];

	//std::cout << "Best order gc:\t" << best_gc << std::endl;
	return { best_order, best_gc };




}


