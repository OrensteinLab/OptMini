#include "dfs.h"

void dfs_initial_wrapper(uint64_t kmer, std::unordered_set<uint64_t> kmers, uint32_t lev, uint32_t w, uint32_t k, std::vector<uint64_t>& gc) {
	uint64_t kmer_mask = (static_cast<uint64_t>(1) << k) - 1;
	uint32_t levels_left = w - lev;
	dfs_initial(kmer, kmers, levels_left, kmer_mask, gc);
}


void dfs_initial(uint64_t kmer, std::unordered_set<uint64_t> kmers, uint32_t levels_left, uint64_t kmer_mask, std::vector<uint64_t>& gc) {
	bool flag = kmers.find(kmer) == kmers.end();
	if (levels_left == 0) {
		if (flag) {
			gc[kmer] += 1;
		}
	}
	else {
		if (flag) {
			kmers.insert(kmer);
		}
		uint64_t next_kmer = (kmer << 1) & kmer_mask;
		dfs_initial(next_kmer, kmers, levels_left - 1, kmer_mask, gc);
		dfs_initial(next_kmer | 1, kmers, levels_left - 1, kmer_mask, gc);
		if (flag) {
			kmers.erase(kmer);
		}
	}
}


uint64_t dfs_prefix_wrapper(uint64_t kmer, uint64_t mask, uint32_t lev, uint32_t w, uint32_t k) {
	uint64_t kmer_mask = (static_cast<uint64_t>(1) << k) - 1;
	uint32_t levels_left = w - lev;
	return dfs_prefix(kmer, mask, levels_left, kmer_mask);

}




uint64_t dfs_prefix(uint64_t kmer, uint64_t mask, uint32_t levels_left, uint64_t kmer_mask) {
	// in case kmer already has a rank
	if ((mask >> kmer) & 1) {
		return 0;
	}
	if (levels_left == 0) {
		return 1;
	}
	uint64_t next_kmer = (kmer << 1) & kmer_mask;
	return dfs_prefix(next_kmer, mask, levels_left - 1, kmer_mask) + dfs_prefix(next_kmer | 1, mask, levels_left - 1, kmer_mask);
}


uint64_t dfs_suffix_wrapper(uint64_t kmer, uint64_t mask, uint32_t lev, uint32_t w, uint32_t k) {
	uint64_t left_bit = static_cast<uint64_t>(1) << (k - 1);
	uint32_t levels_left = w - lev;
	return dfs_suffix(kmer, mask, levels_left, left_bit);
}


uint64_t dfs_suffix(uint64_t kmer, uint64_t mask, uint32_t levels_left, uint64_t left_bit) {
	// in case kmer already has a rank
	if ((mask >> kmer) & 1) {
		return 0;
	}
	if (levels_left == 0) {
		return 1;
	}
	uint64_t previous_kmer = kmer >> 1;
	return dfs_suffix(previous_kmer, mask, levels_left - 1, left_bit) + dfs_suffix(previous_kmer | left_bit, mask, levels_left - 1, left_bit);

}






std::pair<uint64_t, uint64_t> dij_exhaust_dfs(uint32_t w, uint32_t k, uint64_t n_charged) {



	std::cout << "w:\t" << w << "\tk:\t" << k << std::endl;



	uint64_t n_kmers = static_cast<uint64_t>(1) << k;
	uint64_t kmer_mask = n_kmers - 1;
	//std::cout << "mask:\t" << kmer_mask << std::endl;

	std::vector<uint64_t> gc(n_kmers, static_cast<uint64_t>(1) << w);
	for (uint64_t prefix = 0; prefix < n_kmers; prefix++)
	{
		uint64_t next_kmer = (prefix << 1) & kmer_mask;
		dfs_initial_wrapper(next_kmer, { prefix }, 1, w, k, gc);
		dfs_initial_wrapper(next_kmer | 1, { prefix }, 1, w, k, gc);
	}

	// #key is a mask of a subset, value is gc + correction for 0..0 and 1..1
	std::unordered_map<uint64_t, uint64_t> mask_to_gc;

	uint64_t n_kmers_starting_with_0 = static_cast<uint64_t>(1) << (k - 1);

	mask_to_gc[1] = gc[0] + 1; // only need correction for 11111
	for (uint64_t prefix = 1; prefix < n_kmers_starting_with_0; prefix++)
	{
		mask_to_gc[static_cast<uint64_t>(1) << prefix] = gc[prefix] + 2;
	}


	// Set to track seen masks and their gc values
	std::unordered_map<uint64_t, uint64_t> seen_masks; // Maps mask to its best gc_value


	// Min-heap to track gc and masks
	std::priority_queue<std::pair<uint64_t, uint64_t>,
		std::vector<std::pair<uint64_t, uint64_t>>,
		MinHeapComparator>
		minHeap;

	// Populate the heap with mask_to_gc
	for (const auto& [mask, gc_value] : mask_to_gc) {
		minHeap.emplace(gc_value, mask);
	}

	// update the seen masks
	for (const auto& [mask, gc_value] : mask_to_gc) {
		seen_masks[mask] = gc_value;
	}


	uint64_t iteration = 0;

	// Process the heap
	while (!minHeap.empty()) {

		if (iteration % 1'000'000 == 0) {
			std::cout << "[Cleanup] iteration " << iteration
				<< " — cleaning duplicates in heap (" << minHeap.size() << " entries)\n";

			// Temporary map to keep only lowest gc_value per mask
			std::unordered_map<uint64_t, uint64_t> unique_masks;
			while (!minHeap.empty()) {
				auto [val, m] = minHeap.top();
				minHeap.pop();
				auto it = unique_masks.find(m);
				if (it == unique_masks.end() || val < it->second)
					unique_masks[m] = val;
			}

			// Rebuild heap and update seen_masks
			for (const auto& [m, val] : unique_masks) {
				minHeap.emplace(val, m);
				seen_masks[m] = val;
			}

			std::cout << "[Cleanup] after dedup — heap_size=" << minHeap.size() << std::endl;
		}




		auto [gc_value, mask] = minHeap.top(); // Get the minimum gc and its mask
		minHeap.pop();
		iteration++;

		// Skip if this gc_value is outdated
		if (seen_masks[mask] != gc_value) {
			continue;
		}



		if (iteration % 1000 == 0) {
			uint64_t heap_size = minHeap.size();
			std::cout << "iteration:\t" << iteration << "\tgc:\t" << gc_value << "\theap_size:\t" << heap_size << std::endl;
		}


		bool can_be_uhs = true;
		for (uint64_t kmer = 0; kmer < n_kmers; kmer++)
		{
			uint64_t temp_value = gc_value;
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
				uint64_t prefix_to_add = dfs_prefix_wrapper(next_kmer, mask, 1, w, k) + dfs_prefix_wrapper(next_kmer | 1, mask, 1, w, k);
				uint64_t suffix_to_add = dfs_suffix_wrapper(previous_kmer, new_mask, 1, w, k) + dfs_suffix_wrapper(previous_kmer | (static_cast<uint64_t>(1) << (k - 1)), new_mask, 1, w, k);
				temp_value += prefix_to_add + suffix_to_add;
				///std::cout << "temp_value:\t" << temp_value << "\tn_charged:\t" << n_charged << std::endl;
				///std::cout << "prefix_to_add:\t" << prefix_to_add << "\tsuffix_to_add:\t" << suffix_to_add << std::endl;

				// TODO: can be >= too 
				if (temp_value > n_charged) {
					can_be_uhs = false;
					continue;
				}


				if (prefix_to_add > 0 || suffix_to_add > 0) {
					can_be_uhs = false;
				}


				// Check if we have seen this mask before
				auto it = seen_masks.find(new_mask);
				if (it != seen_masks.end()) {
					if (temp_value < it->second) {
						// Update the gc_value in the map
						seen_masks[new_mask] = temp_value;
						// Add to the heap
						minHeap.emplace(temp_value, new_mask);
					}
					// If not lower, do nothing (skip adding to heap)
				}
				else {
					// Add to the heap and seen_masks
					minHeap.emplace(temp_value, new_mask);
					seen_masks[new_mask] = temp_value;
				}

			}
		}
		//std::cout << "Can be UHS:\t" << can_be_uhs << std::endl;
		if (can_be_uhs) {
			// check if the mask is UHS
			if (is_uhs_dfs_wrapper(mask, w, k)) {
				std::cout << "UHS found" << std::endl;
				std::cout << "mask:\t" << mask << "\tgc:\t" << gc_value << std::endl;
				return std::make_pair(mask, gc_value);
			}
		}

	}
	return std::make_pair(0, 0);


}





// ===================== FINAL SINGLE-THREAD POSTCHECK =====================
static void process_remaining_heap_single_dfs(
	std::vector<std::pair<uint64_t, uint64_t>>& globalHeap,
	std::array<SeenShard, N_SHARDS>& seen_shards,
	uint32_t w, uint32_t k, uint64_t n_charged,
	std::atomic<uint64_t>& best_gc_found,
	std::atomic<uint64_t>& best_mask_found)
{
	const uint64_t n_kmers = 1ULL << k;
	const uint64_t kmer_mask = n_kmers - 1ULL;
	size_t processed = 0;

	while (!globalHeap.empty()) {
		std::pop_heap(globalHeap.begin(), globalHeap.end(), MinHeapComparator{});
		auto [gc_value, mask] = globalHeap.back();
		globalHeap.pop_back();
		++processed;

		if (gc_value >= best_gc_found.load())
			break; // nothing better can appear

		bool can_be_uhs = true;
		for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
			if ((mask >> kmer) & 1ULL) continue;

			uint64_t tmp = gc_value;
			if (kmer == 0 || kmer == (n_kmers - 1)) tmp -= 1;

			uint64_t new_mask = mask | (1ULL << kmer);
			uint64_t next_kmer = (kmer << 1) & kmer_mask;
			uint64_t prev_kmer = kmer >> 1;

			uint64_t prefix_to_add =
				dfs_prefix_wrapper(next_kmer, mask, 1, w, k) +
				dfs_prefix_wrapper(next_kmer | 1ULL, mask, 1, w, k);
			uint64_t suffix_to_add =
				dfs_suffix_wrapper(prev_kmer, new_mask, 1, w, k) +
				dfs_suffix_wrapper(prev_kmer | (1ULL << (k - 1)), new_mask, 1, w, k);

			tmp += prefix_to_add + suffix_to_add;
			if (tmp > n_charged) { can_be_uhs = false; continue; }
			if (prefix_to_add > 0 || suffix_to_add > 0) can_be_uhs = false;
		}

		if (can_be_uhs && is_uhs_dfs_wrapper(mask, w, k)) {
			if (gc_value < best_gc_found.load()) {
				best_gc_found.store(gc_value);
				best_mask_found.store(mask);
				std::cout << "[Post-check] Better UHS found: gc=" << gc_value
					<< " mask=" << mask << std::endl;
			}
		}
	}

	std::cout << "[Post-check] scanned " << processed
		<< " items, final best_gc=" << best_gc_found.load() << std::endl;
}

// ===================== MAIN =====================
std::pair<uint64_t, uint64_t>
dij_exhaust_parallel_shared_heap(uint32_t w, uint32_t k, uint64_t n_charged, int N)
{
	if (N < 1) N = 1;
	const uint64_t n_kmers = 1ULL << k;
	const uint64_t kmer_mask = n_kmers - 1ULL;

	// ---- Initial GC via DFS ----
	std::vector<uint64_t> gc(n_kmers, 1ULL << w);
	for (uint64_t prefix = 0; prefix < n_kmers; ++prefix) {
		uint64_t next_kmer = (prefix << 1) & kmer_mask;
		dfs_initial_wrapper(next_kmer, { prefix }, 1, w, k, gc);
		dfs_initial_wrapper(next_kmer | 1ULL, { prefix }, 1, w, k, gc);
	}

	std::array<SeenShard, N_SHARDS> seen_shards;
	for (auto& s : seen_shards) s.map.reserve(200'000);

	std::vector<std::pair<uint64_t, uint64_t>> globalHeap;
	globalHeap.reserve(1'000'000);
	std::mutex heap_mutex;
	std::condition_variable heap_cv;

	std::atomic<bool> stop_requested{ false };
	std::atomic<uint64_t> best_gc_found{ std::numeric_limits<uint64_t>::max() };
	std::atomic<uint64_t> best_mask_found{ 0 };
	std::atomic<bool> uhs_found{ false };
	std::atomic<int>  active_processing{ 0 };
	std::atomic<uint64_t> total_popped{ 0 };

	// ---- Seed ----
	{
		const uint64_t n_kmers_starting_with_0 = 1ULL << (k - 1);
		seen_set(seen_shards, 1ULL, gc[0] + 1);
		globalHeap.emplace_back(gc[0] + 1, 1ULL);
		for (uint64_t prefix = 1; prefix < n_kmers_starting_with_0; ++prefix) {
			uint64_t m = (1ULL << prefix);
			uint64_t v = gc[prefix] + 2;
			seen_set(seen_shards, m, v);
			globalHeap.emplace_back(v, m);
		}
		std::make_heap(globalHeap.begin(), globalHeap.end(), MinHeapComparator{});
	}

	// ---- Worker ----
	auto worker_fn = [&](int worker_id)
		{
			std::vector<std::pair<uint64_t, uint64_t>> local_push;
			local_push.reserve(MAX_LOCAL_PUSH_BATCH);
			std::vector<std::pair<uint64_t, uint64_t>> batch;
			batch.reserve(WORKLOAD_PER_WORKER);

			while (true) {
				if (stop_requested.load()) break;

				{
					std::unique_lock<std::mutex> lk(heap_mutex);
					heap_cv.wait(lk, [&] { return stop_requested.load() || !globalHeap.empty(); });
					if (stop_requested.load()) break;

					size_t count = 0;
					while (!globalHeap.empty() && count < WORKLOAD_PER_WORKER) {
						std::pop_heap(globalHeap.begin(), globalHeap.end(), MinHeapComparator{});
						batch.push_back(globalHeap.back());
						globalHeap.pop_back();
						++count;
					}
				}

				if (batch.empty()) continue;
				active_processing.fetch_add(1);
				total_popped.fetch_add(batch.size());

				for (auto& item : batch) {
					const uint64_t gc_value = item.first;
					const uint64_t mask = item.second;

					uint64_t cur_seen = 0;
					if (seen_lookup(seen_shards, mask, cur_seen) && cur_seen != gc_value)
						continue;

					bool can_be_uhs = true;
					for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
						if ((mask >> kmer) & 1ULL) continue;

						uint64_t tmp = gc_value;
						if (kmer == 0 || kmer == (n_kmers - 1)) tmp -= 1;

						const uint64_t new_mask = mask | (1ULL << kmer);
						const uint64_t next_kmer = (kmer << 1) & kmer_mask;
						const uint64_t prev_kmer = kmer >> 1;

						const uint64_t prefix_to_add =
							dfs_prefix_wrapper(next_kmer, mask, 1, w, k) +
							dfs_prefix_wrapper(next_kmer | 1ULL, mask, 1, w, k);
						const uint64_t suffix_to_add =
							dfs_suffix_wrapper(prev_kmer, new_mask, 1, w, k) +
							dfs_suffix_wrapper(prev_kmer | (1ULL << (k - 1)), new_mask, 1, w, k);

						tmp += prefix_to_add + suffix_to_add;
						if (tmp > n_charged) { can_be_uhs = false; continue; }
						if (prefix_to_add > 0 || suffix_to_add > 0) can_be_uhs = false;

						if (seen_update_if_better(seen_shards, new_mask, tmp)) {
							local_push.emplace_back(tmp, new_mask);
							if (local_push.size() >= MAX_LOCAL_PUSH_BATCH) {
								std::lock_guard<std::mutex> lk(heap_mutex);
								for (auto& p : local_push) {
									globalHeap.push_back(p);
									std::push_heap(globalHeap.begin(), globalHeap.end(), MinHeapComparator{});
								}
								local_push.clear();
								heap_cv.notify_all();
							}
						}
					}

					if (can_be_uhs && is_uhs_dfs_wrapper(mask, w, k)) {
						uint64_t old_best = best_gc_found.load();
						if (gc_value < old_best) {
							best_gc_found.store(gc_value);
							best_mask_found.store(mask);
						}
						uhs_found.store(true);
						stop_requested.store(true);
						heap_cv.notify_all();
					}
				}

				if (!local_push.empty()) {
					std::lock_guard<std::mutex> lk(heap_mutex);
					for (auto& p : local_push) {
						globalHeap.push_back(p);
						std::push_heap(globalHeap.begin(), globalHeap.end(), MinHeapComparator{});
					}
					local_push.clear();
					heap_cv.notify_all();
				}

				batch.clear();
				active_processing.fetch_sub(1);
			}
		};

	std::vector<std::thread> workers;
	for (int i = 0; i < N; ++i)
		workers.emplace_back(worker_fn, i);

	// ---- Main cleanup loop ----
	while (true) {
		std::this_thread::sleep_for(std::chrono::milliseconds(CLEANUP_PERIOD_MS));

		{
			std::lock_guard<std::mutex> lk(heap_mutex);
			rebuild_minheap_dedup(globalHeap, seen_shards);
			std::cout << "[Cleanup] after dedup — heap_size=" << globalHeap.size() << std::endl;
		}
		heap_cv.notify_all();

		if (VERBOSE_STATUS) {
			std::lock_guard<std::mutex> lk(heap_mutex);
			std::cout << "[Status] popped=" << total_popped.load()
				<< " heap=" << globalHeap.size()
				<< " active=" << active_processing.load()
				<< " stop=" << stop_requested.load()
				<< " min heap gc=" << (globalHeap.empty() ? 0 : globalHeap.front().first) << std::endl;
		}

		bool done = false;
		{
			std::lock_guard<std::mutex> lk(heap_mutex);
			const bool heap_empty = globalHeap.empty();
			const bool active0 = (active_processing.load() == 0);
			if ((stop_requested.load() && active0) || (heap_empty && active0))
				done = true;
		}
		if (done) break;
	}

	stop_requested.store(true);
	heap_cv.notify_all();
	for (auto& t : workers) t.join();

	// ---- Single-thread post-check ----
	if (uhs_found.load()) {
		std::cout << "[Post-check] continuing with one thread up to GC ≤ "
			<< best_gc_found.load() << std::endl;
		process_remaining_heap_single_dfs(globalHeap, seen_shards, w, k, n_charged,
			best_gc_found, best_mask_found);
	}

	if (best_mask_found.load() != 0) {
		std::cout << "Final UHS: mask=" << best_mask_found.load()
			<< " gc=" << best_gc_found.load() << std::endl;
		return { best_mask_found.load(), best_gc_found.load() };
	}
	else {
		std::cout << "No UHS found." << std::endl;
		return { 0,0 };
	}
}




std::vector<uint64_t> get_best_order(uint64_t original_mask, uint64_t best_gc, uint32_t w, uint32_t k) {


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

	uint64_t n_kmers_starting_with_0 = static_cast<uint64_t>(1) << (k - 1);

	if ((original_mask >> 0) & 1) {
		old_mask_to_gc[1] = gc[0] + 1; // only need correction for 11111
		mask_to_best_order[1] = { 0 };
	}

	for (uint64_t prefix = 1; prefix < n_kmers_starting_with_0; prefix++)
	{
		if ((original_mask >> prefix) & 1) {
			old_mask_to_gc[static_cast<uint64_t>(1) << prefix] = gc[prefix] + 2;
			mask_to_best_order[static_cast<uint64_t>(1) << prefix] = { prefix };
		}

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

			if (is_uhs_dfs_wrapper(mask, w, k)) {
				std::cout << "UHS found" << std::endl;
				std::cout << "mask:\t" << mask << "\tgc:\t" << value << std::endl;
				std::cout << "best order:\t";
				for (auto& kmer : mask_to_best_order[mask]) {
					std::cout << kmer << " ";
				}
				std::cout << std::endl;
				return mask_to_best_order[mask];
			}

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
					uint64_t next_kmer = (kmer << 1) & kmer_mask;
					uint64_t previous_kmer = kmer >> 1;
					temp_value += dfs_prefix_wrapper(next_kmer, mask, 1, w, k) + dfs_prefix_wrapper(next_kmer | 1, mask, 1, w, k);
					temp_value += dfs_suffix_wrapper(previous_kmer, new_mask, 1, w, k) + dfs_suffix_wrapper(previous_kmer | (static_cast<uint64_t>(1) << (k - 1)), new_mask, 1, w, k);
					if (temp_value <= best_gc) {
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
			// do an error
			std::cout << "Couldn't reconstruct an order from mask and gc";
			return {};
		}

	}





}



