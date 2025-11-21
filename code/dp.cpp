#include "dp.h"

// Compute GC array using DP (not DFS)
static std::vector<uint64_t> compute_gc_dp(uint32_t w, uint32_t k)
{
	uint64_t n_kmers = 1ULL << k;
	uint64_t kmer_mask = n_kmers - 1ULL;
	std::vector<uint64_t> gc(n_kmers, 1ULL << w);

	for (uint64_t suffix = 0; suffix < n_kmers; ++suffix) {
		std::vector<uint64_t> old_arr(n_kmers, 1ULL);
		std::vector<uint64_t> new_arr(n_kmers, 0ULL);
		for (uint64_t i = 0; i < w; ++i) {
			old_arr[suffix] = 0ULL;
			for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
				uint64_t next_kmer = (kmer << 1) & kmer_mask;
				new_arr[next_kmer] += old_arr[kmer];
				new_arr[next_kmer | 1] += old_arr[kmer];
			}
			std::fill(old_arr.begin(), old_arr.end(), 0ULL);
			std::swap(old_arr, new_arr);
			std::fill(new_arr.begin(), new_arr.end(), 0ULL);
		}
		gc[suffix] += old_arr[suffix];
	}
	return gc;
}

// Prefix charged via DP
static uint64_t prefix_charged_dp(uint64_t mask, uint64_t kmer_to_add,
	uint32_t w, uint32_t k)
{
	uint64_t n_kmers = 1ULL << k;
	uint64_t kmer_mask = n_kmers - 1ULL;
	std::vector<uint64_t> old_arr(n_kmers, 0ULL);
	std::vector<uint64_t> new_arr(n_kmers, 0ULL);
	old_arr[kmer_to_add] = 1ULL;

	for (uint64_t i = 0; i < w; ++i) {
		for (uint64_t curr_kmer = 0; curr_kmer < n_kmers; ++curr_kmer) {
			if (!((mask >> curr_kmer) & 1ULL)) {
				uint64_t prev_kmer = curr_kmer >> 1;
				new_arr[curr_kmer] = old_arr[prev_kmer] +
					old_arr[prev_kmer | (1ULL << (k - 1))];
			}
		}
		std::swap(old_arr, new_arr);
		std::fill(new_arr.begin(), new_arr.end(), 0ULL);
	}
	return std::accumulate(old_arr.begin(), old_arr.end(), 0ULL);
}

// Suffix charged via DP
static uint64_t suffix_charged_dp(uint64_t new_mask, uint64_t kmer_to_add,
	uint32_t w, uint32_t k)
{
	uint64_t n_kmers = 1ULL << k;
	uint64_t kmer_mask = n_kmers - 1ULL;
	std::vector<uint64_t> old_arr(n_kmers, 0ULL);
	std::vector<uint64_t> new_arr(n_kmers, 0ULL);
	old_arr[kmer_to_add] = 1ULL;

	for (uint64_t i = 0; i < w; ++i) {
		for (uint64_t curr_kmer = 0; curr_kmer < n_kmers; ++curr_kmer) {
			if (!((new_mask >> curr_kmer) & 1ULL)) {
				uint64_t next_kmer = (curr_kmer << 1) & kmer_mask;
				new_arr[curr_kmer] = old_arr[next_kmer] +
					old_arr[next_kmer | 1ULL];
			}
		}
		std::swap(old_arr, new_arr);
		std::fill(new_arr.begin(), new_arr.end(), 0ULL);
	}
	return std::accumulate(old_arr.begin(), old_arr.end(), 0ULL);
}




std::pair<uint64_t, uint64_t> dij_exhaust_dp(uint32_t w, uint32_t k, uint64_t n_charged) {



	std::cout << "w:\t" << w << "\tk:\t" << k << std::endl;



	uint64_t n_kmers = static_cast<uint64_t>(1) << k;
	uint64_t kmer_mask = n_kmers - 1;
	//std::cout << "mask:\t" << kmer_mask << std::endl;

	// initialized with prefix GCs
	std::vector<uint64_t> gc(n_kmers, static_cast<uint64_t>(1) << w);

	for (uint64_t suffix = 0; suffix < n_kmers; suffix++)
	{
		// create a uint64_t array of size 2^k filled with ones called old
		std::vector<uint64_t> old_arr(n_kmers, static_cast<uint64_t>(1));
		std::vector<uint64_t> new_arr(n_kmers, static_cast<uint64_t>(0));
		for (uint64_t i = 0; i < w; i++)
		{
			old_arr[suffix] = 0; // discard windows that arent suffix charged
			for (uint64_t kmer = 0; kmer < n_kmers; kmer++)
			{
				uint64_t next_kmer = (kmer << 1) & kmer_mask;
				new_arr[next_kmer] += old_arr[kmer];
				new_arr[next_kmer | 1] += old_arr[kmer];
			}
			// reset old_arr to zeros
			old_arr = std::vector<uint64_t>(n_kmers, static_cast<uint64_t>(0));
			// swap old and new
			std::swap(old_arr, new_arr);
		}
		// add to gc
		gc[suffix] += old_arr[suffix];

		//std::cout << "suffix:\t" << std::bitset<64>(suffix) << "\tgc:\t" << gc[suffix] << std::endl;

	}


	// #key is a mask of a subset, value is gc + correction for 0..0 and 1..1
	std::unordered_map<uint64_t, uint64_t> mask_to_gc;

	uint64_t n_kmers_starting_with_0 = static_cast<uint64_t>(1) << (k - 1);



	// only need correction for 11111 for the 00000 window ( since usually need 2 corrections)
	// we use 1 because its basically  1<<0 (the mask containing the 00000 kmer)
	mask_to_gc[1] = gc[0] + 1;
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
		for (uint64_t kmer_to_be_added_to_mask = 0; kmer_to_be_added_to_mask < n_kmers; kmer_to_be_added_to_mask++)
		{
			uint64_t temp_value = gc_value;
			// if kmer isn't in the mask
			if (!((mask >> kmer_to_be_added_to_mask) & 1))
			{
				// handle 1111 or 00000
				if (kmer_to_be_added_to_mask == 0 || kmer_to_be_added_to_mask == (n_kmers - 1))
				{
					temp_value -= 1;
				}

				uint64_t new_mask = mask | (static_cast<uint64_t>(1) << kmer_to_be_added_to_mask);


				//std::cout << "Old mask:\t" << std::bitset<64>(mask) << "\tKmer to add:\t" << kmer_to_be_added_to_mask << "\tNew mask:\t" << std::bitset<64>(new_mask) << std::endl;




				// calc prefix charged//////////////////////////////////////////////////////////
				std::vector<uint64_t> old_arr2(n_kmers, static_cast<uint64_t>(0));
				std::vector<uint64_t> new_arr2(n_kmers, static_cast<uint64_t>(0));
				old_arr2[kmer_to_be_added_to_mask] = 1;

				for (uint64_t i = 0; i < w; i++)
				{
					for (uint64_t curr_kmer = 0; curr_kmer < n_kmers; curr_kmer++)
					{
						// only proceed if curr_kmer is not in the mask (otherwise it wont be prefix charged)
						if (!((mask >> curr_kmer) & 1))
						{
							uint64_t prev_kmer = (curr_kmer >> 1);
							new_arr2[curr_kmer] = (old_arr2[prev_kmer] + old_arr2[prev_kmer | (static_cast<uint64_t>(1) << (k - 1))]);
						}
					}
					// reset old_arr to zeros
					std::swap(old_arr2, new_arr2);
					std::fill(new_arr2.begin(), new_arr2.end(), 0ULL);

				}

				// old_arr2 sum should be added to temp_value

				uint64_t prefix_to_add = std::accumulate(old_arr2.begin(), old_arr2.end(), 0ULL);
				temp_value += prefix_to_add;



				// calc suffix charged///////////////////////////////////////////////////////////////////////////////
				std::fill(new_arr2.begin(), new_arr2.end(), 0ULL);
				std::fill(old_arr2.begin(), old_arr2.end(), 0ULL);
				old_arr2[kmer_to_be_added_to_mask] = 1;

				for (uint64_t i = 0; i < w; i++)
				{
					for (uint64_t curr_kmer = 0; curr_kmer < n_kmers; curr_kmer++)
					{
						// only proceed if curr_kmer is not in the new mask (otherwise it wont be suffix charged)
						if (!((new_mask >> curr_kmer) & 1))
						{
							uint64_t next_kmer = (curr_kmer << 1) & kmer_mask;
							new_arr2[curr_kmer] = (old_arr2[next_kmer] + old_arr2[next_kmer | 1]);
						}
					}
					std::swap(old_arr2, new_arr2);
					std::fill(new_arr2.begin(), new_arr2.end(), 0ULL);

				}

				uint64_t suffix_to_add = std::accumulate(old_arr2.begin(), old_arr2.end(), 0ULL);
				temp_value += suffix_to_add;






				if (prefix_to_add > 0 || suffix_to_add > 0) {
					can_be_uhs = false;
				}



				if (temp_value > n_charged) {
					continue;
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
			// print mask
			//std::cout << "Mask : \t" << std::bitset<64>(mask) << std::endl;

			// check if the mask is UHS
			if (is_uhs_dp_wrapper(mask, w, k)) {
				std::cout << "UHS found" << std::endl;
				std::cout << "mask:\t" << mask << "\tgc:\t" << gc_value << std::endl;
				return std::make_pair(mask, gc_value);
			}
		}

	}
	return std::make_pair(0, 0);






}




std::pair<uint64_t, uint64_t>
dij_exhaust_parallel_shared_heap_dp(uint32_t w, uint32_t k, uint64_t n_charged, int N)
{
	if (N < 1) N = 1;
	const uint64_t n_kmers = 1ULL << k;
	const uint64_t kmer_mask = n_kmers - 1ULL;

	std::vector<uint64_t> gc = compute_gc_dp(w, k);
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
						const uint64_t prefix_to_add = prefix_charged_dp(mask, kmer, w, k);
						const uint64_t suffix_to_add = suffix_charged_dp(new_mask, kmer, w, k);

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

					if (can_be_uhs && is_uhs_dp_wrapper(mask, w, k)) {
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
				<< " stop=" << stop_requested.load() << std::endl
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
		process_remaining_heap_single_dp(globalHeap, seen_shards, w, k, n_charged,
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




// same but dont add to mask if it doesnt increase charged contexts AND doesnt increase unchraged contexts (useless kmer)
std::pair<uint64_t, uint64_t>
dij_exhaust_parallel_shared_heap_dp_restrict_masks(uint32_t w, uint32_t k, uint64_t n_charged, int N)
{
	if (N < 1) N = 1;
	const uint64_t n_kmers = 1ULL << k;
	const uint64_t kmer_mask = n_kmers - 1ULL;

	std::vector<uint64_t> gc = compute_gc_dp(w, k);
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
						const uint64_t prefix_to_add = prefix_charged_dp(mask, kmer, w, k);
						const uint64_t suffix_to_add = suffix_charged_dp(new_mask, kmer, w, k);

						// skip k-mers that dont help us cover more contexts
						if (prefix_to_add == 0 && suffix_to_add == 0 && is_kmer_useless_dp(w, k, mask, kmer))
						{
							continue;
						}

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

					if (can_be_uhs && is_uhs_dp_wrapper(mask, w, k)) {
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
				<< " stop=" << stop_requested.load() << std::endl
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
		process_remaining_heap_single_dp_restrict_masks(globalHeap, seen_shards, w, k, n_charged,
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



static void process_remaining_heap_single_dp(
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

			uint64_t prefix_to_add = prefix_charged_dp(mask, kmer, w, k);
			uint64_t suffix_to_add = suffix_charged_dp(new_mask, kmer, w, k);

			tmp += prefix_to_add + suffix_to_add;
			if (tmp > n_charged) { can_be_uhs = false; continue; }
			if (prefix_to_add > 0 || suffix_to_add > 0) can_be_uhs = false;
		}

		if (can_be_uhs && is_uhs_dp_wrapper(mask, w, k)) {
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




static void process_remaining_heap_single_dp_restrict_masks(
	std::vector<std::pair<uint64_t, uint64_t>>& globalHeap,
	std::array<SeenShard, N_SHARDS>& seen_shards,
	uint32_t w, uint32_t k, uint64_t n_charged,
	std::atomic<uint64_t>& best_gc_found,
	std::atomic<uint64_t>& best_mask_found)
{
	const uint64_t n_kmers = 1ULL << k;
	const uint64_t kmer_mask = n_kmers - 1ULL;
	size_t processed = 0;

	// rebuild heap just in case
	std::make_heap(globalHeap.begin(), globalHeap.end(), MinHeapComparator{});

	while (!globalHeap.empty()) {
		std::pop_heap(globalHeap.begin(), globalHeap.end(), MinHeapComparator{});
		auto [gc_value, mask] = globalHeap.back();
		globalHeap.pop_back();
		++processed;

		// --- Stop when GC ≥ best found ---
		if (gc_value >= best_gc_found.load())
			break;

		uint64_t cur_seen = 0;
		if (seen_lookup(seen_shards, mask, cur_seen) && cur_seen != gc_value)
			continue;

		bool can_be_uhs = true;

		for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
			if ((mask >> kmer) & 1ULL) continue; // already in mask

			uint64_t tmp = gc_value;
			if (kmer == 0 || kmer == (n_kmers - 1)) tmp -= 1ULL;

			const uint64_t new_mask = mask | (1ULL << kmer);
			const uint64_t prefix_to_add = prefix_charged_dp(mask, kmer, w, k);
			const uint64_t suffix_to_add = suffix_charged_dp(new_mask, kmer, w, k);

			// --- Skip useless k-mers ---
			if (prefix_to_add == 0 && suffix_to_add == 0 &&
				is_kmer_useless_dp(w, k, mask, kmer))
				continue;

			tmp += prefix_to_add + suffix_to_add;
			if (tmp > n_charged) { can_be_uhs = false; continue; }
			if (prefix_to_add > 0 || suffix_to_add > 0) can_be_uhs = false;

			// --- Update seen and push if valuable ---
			if (seen_update_if_better(seen_shards, new_mask, tmp)) {
				if (tmp < best_gc_found.load()) {
					globalHeap.emplace_back(tmp, new_mask);
					std::push_heap(globalHeap.begin(), globalHeap.end(), MinHeapComparator{});
				}
			}
		}

		// --- Check UHS condition ---
		if (can_be_uhs && is_uhs_dp_wrapper(mask, w, k)) {
			uint64_t old_best = best_gc_found.load();
			if (gc_value < old_best) {
				best_gc_found.store(gc_value);
				best_mask_found.store(mask);
				std::cout << "[Post-check] Better UHS found: gc=" << gc_value
					<< " mask=" << mask << std::endl;

				// prune entries >= new best
				globalHeap.erase(
					std::remove_if(globalHeap.begin(), globalHeap.end(),
						[&](const auto& p) { return p.first >= gc_value; }),
					globalHeap.end());
				std::make_heap(globalHeap.begin(), globalHeap.end(), MinHeapComparator{});
			}
		}
	}

	std::cout << "[Post-check] scanned " << processed
		<< " items, final best_gc=" << best_gc_found.load()
		<< " heap_left=" << globalHeap.size() << std::endl;
}










std::vector<uint64_t> get_best_order_dp(uint64_t original_mask, uint64_t best_gc, uint32_t w, uint32_t k)
{
	const uint64_t n_kmers = 1ULL << k;
	const uint64_t kmer_mask = n_kmers - 1ULL;

	// ---------- Compute base GC table via DP ----------
	std::vector<uint64_t> gc(n_kmers, 1ULL << w);
	for (uint64_t suffix = 0; suffix < n_kmers; ++suffix) {
		std::vector<uint64_t> old_arr(n_kmers, 1ULL);
		std::vector<uint64_t> new_arr(n_kmers, 0ULL);

		for (uint64_t i = 0; i < w; ++i) {
			old_arr[suffix] = 0; // discard windows that aren't suffix charged
			for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
				uint64_t next_kmer = (kmer << 1) & kmer_mask;
				new_arr[next_kmer] += old_arr[kmer];
				new_arr[next_kmer | 1] += old_arr[kmer];
			}
			std::swap(old_arr, new_arr);
			std::fill(new_arr.begin(), new_arr.end(), 0ULL);
		}
		gc[suffix] += old_arr[suffix];
	}

	// ---------- Initialize with singleton masks ----------
	std::unordered_map<uint64_t, uint64_t> old_mask_to_gc;
	std::unordered_map<uint64_t, std::vector<uint64_t>> mask_to_best_order;

	const uint64_t n_kmers_starting_with_0 = 1ULL << (k - 1);

	if ((original_mask >> 0) & 1ULL) {
		old_mask_to_gc[1ULL] = gc[0] + 1;
		mask_to_best_order[1ULL] = { 0 };
	}

	for (uint64_t prefix = 1; prefix < n_kmers_starting_with_0; ++prefix) {
		if ((original_mask >> prefix) & 1ULL) {
			uint64_t m = 1ULL << prefix;
			old_mask_to_gc[m] = gc[prefix] + 2;
			mask_to_best_order[m] = { prefix };
		}
	}

	// ---------- Main expansion cycle ----------
	for (uint64_t rank = 0; rank <= n_kmers; ++rank) {
		std::unordered_map<uint64_t, uint64_t> mask_to_gc;

		for (auto& [mask, value] : old_mask_to_gc) {
			if (is_uhs_dp_wrapper(mask, w, k)) {
				std::cout << "UHS found\n";
				std::cout << "mask:\t" << mask << "\tgc:\t" << value << "\n";
				std::cout << "best order:\t";
				for (auto& kmer : mask_to_best_order[mask])
					std::cout << kmer << " ";
				std::cout << std::endl;
				return mask_to_best_order[mask];
			}

			for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
				if (!((original_mask >> kmer) & 1ULL))
					continue;

				uint64_t temp_value = value;
				if (!((mask >> kmer) & 1ULL)) {
					if (kmer == 0 || kmer == (n_kmers - 1))
						temp_value -= 1;

					uint64_t new_mask = mask | (1ULL << kmer);

					// -------- prefix charged DP ------------
					std::vector<uint64_t> old_arr2(n_kmers, 0ULL);
					std::vector<uint64_t> new_arr2(n_kmers, 0ULL);
					old_arr2[kmer] = 1;

					for (uint64_t i = 0; i < w; ++i) {
						for (uint64_t curr_kmer = 0; curr_kmer < n_kmers; ++curr_kmer) {
							if (!((mask >> curr_kmer) & 1ULL)) {
								uint64_t prev_kmer = curr_kmer >> 1;
								new_arr2[curr_kmer] =
									old_arr2[prev_kmer] +
									old_arr2[prev_kmer | (1ULL << (k - 1))];
							}
						}
						std::swap(old_arr2, new_arr2);
						std::fill(new_arr2.begin(), new_arr2.end(), 0ULL);
					}
					uint64_t prefix_to_add = std::accumulate(old_arr2.begin(), old_arr2.end(), 0ULL);
					temp_value += prefix_to_add;

					// -------- suffix charged DP ------------
					std::fill(old_arr2.begin(), old_arr2.end(), 0ULL);
					std::fill(new_arr2.begin(), new_arr2.end(), 0ULL);
					old_arr2[kmer] = 1;

					for (uint64_t i = 0; i < w; ++i) {
						for (uint64_t curr_kmer = 0; curr_kmer < n_kmers; ++curr_kmer) {
							if (!((new_mask >> curr_kmer) & 1ULL)) {
								uint64_t next_kmer = (curr_kmer << 1) & kmer_mask;
								new_arr2[curr_kmer] =
									old_arr2[next_kmer] +
									old_arr2[next_kmer | 1ULL];
							}
						}
						std::swap(old_arr2, new_arr2);
						std::fill(new_arr2.begin(), new_arr2.end(), 0ULL);
					}
					uint64_t suffix_to_add = std::accumulate(old_arr2.begin(), old_arr2.end(), 0ULL);
					temp_value += suffix_to_add;

					// -------- filtering and updating ----------
					if (temp_value <= best_gc) {
						if (mask_to_gc.find(new_mask) == mask_to_gc.end() ||
							temp_value < mask_to_gc[new_mask]) {
							mask_to_gc[new_mask] = temp_value;
							mask_to_best_order[new_mask] = mask_to_best_order[mask];
							mask_to_best_order[new_mask].push_back(kmer);
						}
					}
				}
			}
		}

		old_mask_to_gc = mask_to_gc;

		if (old_mask_to_gc.empty()) {
			std::cout << "Couldn't reconstruct an order from mask and gc\n";
			return {};
		}

		uint64_t max_gc = std::max_element(
			old_mask_to_gc.begin(), old_mask_to_gc.end(),
			[](const auto& a, const auto& b) { return a.second < b.second; })->second;

		uint64_t min_gc = std::min_element(
			old_mask_to_gc.begin(), old_mask_to_gc.end(),
			[](const auto& a, const auto& b) { return a.second < b.second; })->second;

		std::cout << "mask size\t" << (rank + 2)
			<< "\tn_masks:\t" << old_mask_to_gc.size()
			<< "\tmin_gc:\t" << min_gc
			<< "\tmax_gc:\t" << max_gc << std::endl;
	}

	return {};
}


using boost::multiprecision::uint256_t;



// =====================================================
// =============== DP Helper Functions (256) ============
// =====================================================

static std::vector<uint256_t> compute_gc_dp_256(uint32_t w, uint32_t k)
{

	const size_t n = 1ULL << k;          // up to 64
	const uint64_t kmer_mask = n - 1ULL; // 64-bit mask for bit ops

	std::vector<uint256_t> gc(n, uint256_t(1) << w);

	for (size_t suffix = 0; suffix < n; ++suffix) {
		std::vector<uint256_t> old_arr(n, 1);
		std::vector<uint256_t> new_arr(n, 0);

		for (uint32_t i = 0; i < w; ++i) {
			old_arr[suffix] = 0;  // clear suffix-charged contexts

			for (size_t kmer = 0; kmer < n; ++kmer) {
				size_t next_kmer = (kmer << 1) & kmer_mask;
				new_arr[next_kmer] += old_arr[kmer];
				new_arr[next_kmer | 1] += old_arr[kmer];
			}

			std::swap(old_arr, new_arr);
			std::fill(new_arr.begin(), new_arr.end(), uint256_t(0));
		}

		gc[suffix] += old_arr[suffix];
	}

	return gc;
}

uint256_t prefix_charged_dp_256(uint64_t mask, uint64_t kmer_to_add,
	uint32_t w, uint32_t k)
{
	const uint64_t n = 1ULL << k;  // ≤ 64
	const uint64_t biggest_bit = 1ULL << (k - 1);
	std::vector<uint256_t> old_arr(n, 0), new_arr(n, 0);
	old_arr[kmer_to_add] = 1;

	for (uint32_t i = 0; i < w; ++i) {
		for (uint64_t curr = 0; curr < n; ++curr) {
			// only propagate if k-mer is *not* masked
			if (((mask >> curr) & 1ULL) == 0) {
				const uint64_t prev = curr >> 1;
				const uint64_t alt = prev | biggest_bit;
				new_arr[curr] = old_arr[prev] + old_arr[alt];
			}
		}
		uint256_t sum = std::accumulate(new_arr.begin(), new_arr.end(), uint256_t(0));
		if (sum == 0) {
			return 0;
		}
		std::swap(old_arr, new_arr);
		std::fill(new_arr.begin(), new_arr.end(), uint256_t(0));



	}

	return std::accumulate(old_arr.begin(), old_arr.end(), uint256_t(0));
}

using boost::multiprecision::uint256_t;
uint256_t suffix_charged_dp_256(uint64_t new_mask, uint64_t kmer_to_add,
	uint32_t w, uint32_t k)
{
	const uint64_t n = 1ULL << k;          // ≤ 64
	const uint64_t kmer_mask = n - 1ULL; // wraparound for (kmer << 1)
	std::vector<uint256_t> old_arr(n, 0), new_arr(n, 0);
	old_arr[kmer_to_add] = 1;

	for (uint32_t i = 0; i < w; ++i) {
		for (uint64_t curr = 0; curr < n; ++curr) {
			if (((new_mask >> curr) & 1ULL) == 0) {
				const uint64_t next = (curr << 1) & kmer_mask;
				new_arr[curr] = old_arr[next] + old_arr[next | 1ULL];
			}
		}
		uint256_t sum = std::accumulate(new_arr.begin(), new_arr.end(), uint256_t(0));
		if (sum == 0) {
			return 0;
		}
		std::swap(old_arr, new_arr);
		std::fill(new_arr.begin(), new_arr.end(), uint256_t(0));
	}

	return std::accumulate(old_arr.begin(), old_arr.end(), uint256_t(0));
}

uint256_t suffix_charged_dp_256_v2(uint64_t old_mask, uint64_t kmer_to_add,
	uint32_t w, uint32_t k)
{
	const uint64_t n = 1ULL << k;          // ≤ 64
	const uint64_t kmer_mask = n - 1ULL; // wraparound for (kmer << 1)
	const uint64_t biggest_bit = 1ULL << (k - 1);

	// array means how many (assuming last kmer which will be ahead is the kmer_to_add suffix) contexts are suffix charged
	std::vector<uint256_t> old_arr(n, 1), new_arr(n, 0);
	// all branches starting from it are now suffix charged
	old_arr[kmer_to_add] = 0;
	for (uint64_t kmer = 0; kmer < n; ++kmer) {
		if (((old_mask >> kmer) & 1ULL) == 1) {
			old_arr[kmer] = 0;
		}
	}


	// first iteration is w=2 (w=1 is handled by initialization), we need to reach w, since eventually we add 1 at the end and get w+1 -> need w-1 iterations
	for (uint32_t i = 1; i < w; ++i) {
		for (uint64_t curr = 0; curr < n; ++curr) {
			if (((old_mask >> curr) & 1ULL) == 0) {
				const uint64_t prev = curr >> 1;
				new_arr[curr] = old_arr[prev] + old_arr[prev | biggest_bit];
			}
		}
		if (new_arr[kmer_to_add] == 0) {
			return 0;
		}
		else {
			new_arr[kmer_to_add] = 0;
		}
		std::swap(old_arr, new_arr);
		std::fill(new_arr.begin(), new_arr.end(), uint256_t(0));
	}

	uint64_t prev_to_kmer_to_add_1 = kmer_to_add >> 1;
	uint64_t prev_to_kmer_to_add_2 = prev_to_kmer_to_add_1 | biggest_bit;
	uint256_t sum = 0;
	if (prev_to_kmer_to_add_1 != kmer_to_add)
		sum += old_arr[prev_to_kmer_to_add_1];
	if (prev_to_kmer_to_add_2 != kmer_to_add)
		sum += old_arr[prev_to_kmer_to_add_2];

	return sum;
}




std::vector<uint256_t> trimmed_dfs_all_living(uint256_t curr_context, uint64_t mask, uint64_t curr_kmer, uint32_t k, uint32_t w_left) {
	uint64_t n = 1ULL << k;

	if (w_left == 0) {
		// a vector  of curr_context
		return { curr_context };
	}

	uint64_t next_kmer = (curr_kmer << 1) & (n - 1ULL);
	uint64_t next_kmer_2 = next_kmer | 1ULL;

	// make an empty vector
	std::vector<uint256_t> result;

	// if next_kmer is not in set
	if (((mask >> next_kmer) & 1ULL) == 0) {
		uint256_t next_context = curr_context << 1; // shift left
		auto res1 = trimmed_dfs_all_living(next_context, mask, next_kmer, k, w_left-1);
		result.insert(result.end(), res1.begin(), res1.end());
	}

	if (((mask >> next_kmer_2) & 1ULL) == 0) {
		uint256_t next_context = (curr_context << 1) | 1ULL; // shift left and add 1
		auto res1 = trimmed_dfs_all_living(next_context, mask, next_kmer_2, k, w_left - 1);
		result.insert(result.end(), res1.begin(), res1.end());
	}

	return result;



}

void trimmed_dfs_all_living_helper(uint64_t mask, uint32_t k, uint32_t w) {
	uint64_t n = 1ULL << k;
	std::vector<uint256_t> result;

	for (uint64_t kmer = 0; kmer < n; ++kmer) {
		// only for kmers not in the mask
		if (((mask >> kmer) & 1ULL) == 0) {
			auto res1 = trimmed_dfs_all_living(kmer,mask, kmer, k, w);
			result.insert(result.end(), res1.begin(), res1.end());
		}
	}

	// print all contexts

	std::cout << "All living contexts in binary (w+k bits):\n";
	for (auto& ctx : result) {
		std::string s = "";
		uint256_t temp = ctx;
		for (uint32_t i = 0; i < w + k; ++i) {
			if ((temp & 1ULL) == 1ULL) {
				s = "1" + s;
			}
			else {
				s = "0" + s;
			}
			temp = temp >> 1;

			if (i == w+k-7) {
				s = "-" + s;
			}
		}
		std::cout << s << "\n";
	}
}

void analyze_mask(uint64_t mask, uint64_t useless_kmers_mask, uint32_t k, uint32_t w) {
	uint64_t n = 1ULL << k;

	// Dictionary: kmer → total prefix+suffix charges
	std::unordered_map<uint64_t, uint256_t> totals;
	std::unordered_map<uint64_t, uint256_t> prefixes;
	std::unordered_map<uint64_t, uint256_t> suffixes; 
	// List of k-mers that are not useless
	std::vector<uint64_t> potential_to_add;

	for (uint64_t kmer = 0; kmer < n; ++kmer) {
		if ((mask >> kmer) & 1ULL) continue;  // already selected
		if (is_kmer_useless_dp_256(w, k, mask, kmer)) continue;

		potential_to_add.push_back(kmer);

		uint256_t prefix = prefix_charged_dp_256(mask, kmer, w, k);
		uint256_t suffix = suffix_charged_dp_256_v2(mask, kmer, w, k);
		uint256_t total = prefix + suffix;

		totals[kmer] = total;  // store in dictionary
		prefixes[kmer] = prefix;
		suffixes[kmer] = suffix;
	}

	// if sum of total is less than 300 -> run dfs
	uint256_t sum = 0;
	for (auto& p : totals) {
		sum += p.second;
	}

	if (sum < 1000) {
		trimmed_dfs_all_living_helper(mask, k, w);
	}




	// Pretty-print dictionary
	//std::cout << "kmer → charge totals:\n[ ";
	//for (auto& p : totals) {
	//	std::cout << p.second << " ";
	//}
	//std::cout << "]\n";


	// pretty-print other dicts together
	std::cout << "kmer → prefix / suffix charges:\n[";
	for (auto& kmer : potential_to_add) {
		std::cout << "" << prefixes[kmer]
			<< "-" << suffixes[kmer] << " ";
	}
	std::cout << "]\n";

	return;

	// for each potential to add:
	for (size_t i = 0; i < potential_to_add.size(); ++i) {
		uint64_t kmer = potential_to_add[i];
		uint64_t mask_with_kmer = mask | (1ULL << kmer);
		uint256_t total_1 = totals[kmer];
		// check if it dominates another kmer to add
		for (size_t j = 0; j < potential_to_add.size(); ++j) {
			if (i == j) continue;
			uint64_t kmer2 = potential_to_add[j];
			uint256_t total_2 = totals[kmer2];

			if (total_1 <= total_2) {
				uint256_t prefix = prefix_charged_dp_256(mask_with_kmer, kmer2, w, k);
				uint256_t suffix = suffix_charged_dp_256_v2(mask_with_kmer, kmer2, w, k);
				uint256_t total = prefix + suffix;
				// if picking kmer 1 will always decrease the total of kmer 2 by more-> no reason to ever add kmer 2 before kmer 1
				if (total_1 + total <= total_2 && total_2 < 5) {
					std::cout << "kmer " << kmer << " (total " << total_1 << ") dominates kmer "
						<< kmer2 << " (total " << total_2 << ")\n";
				}
			}


		}
	}




}

// known: kmer charges 1 prefix + 0 suffixes
// want to check: does k-mer appear in exactly 1 charged context
// how?
// if we look at the set of all kmers unpicked, then this kmer either has 1 or 2 kmers following it OR 1 or 2 kmers preceeding it-> never both!
// cycles are at most length k (assume w>k) in the DBJ graph (min cycles), therefore arrows going into the cycles will always need to resolve the cycle anyways, so it makes no sense to add the kmer 
bool is_kmer_appearing_once(uint64_t mask, uint64_t kmer, uint32_t k, uint32_t w) {

	uint64_t n_kmers = 1ULL << k;
	uint64_t next_kmer = (kmer << 1) & (n_kmers - 1ULL);
	uint64_t next_kmer_2 = next_kmer | 1ULL;
	uint64_t prev_kmer = kmer >> 1;
	uint64_t prev_kmer_2 = prev_kmer | (1ULL << (k - 1));

	// check if those kmers are also not picked yet
	bool has_next_1 = ((mask >> next_kmer) & 1ULL) == 0;
	bool has_next_2 = ((mask >> next_kmer_2) & 1ULL) == 0;
	bool has_prev_1 = ((mask >> prev_kmer) & 1ULL) == 0;
	bool has_prev_2 = ((mask >> prev_kmer_2) & 1ULL) == 0;

	bool has_any_next = has_next_1 || has_next_2;
	bool has_any_prev = has_prev_1 || has_prev_2;

	// if doesnt have prev or doesnt have next -> appears once
	if (!has_any_prev || !has_any_next) {
		return true;
	}
	else {
		return false;
	}
}


struct HeapNode256 {
	boost::multiprecision::uint256_t gc;
	uint64_t mask;
	uint64_t useless_mask;
};
// =====================================================
// =================== Dij Exhaust (256) ================
// =====================================================
std::pair<uint64_t, boost::multiprecision::uint256_t>
dij_exhaust_dp_256_time_limit(uint32_t w, uint32_t k,
	boost::multiprecision::uint256_t n_charged,
	double time_limit_minutes, uint32_t clean_every_n)
{
	using boost::multiprecision::uint256_t;
	using clock = std::chrono::steady_clock;

	const auto start_time = clock::now();
	const auto time_limit = std::chrono::duration<double>(time_limit_minutes * 60.0);

	std::cout << "[timed] w:\t" << w << "\tk:\t" << k
		<< "\tlimit:\t" << time_limit_minutes << " min" << std::endl;

	const uint64_t n = 1ULL << k;  // k ≤ 6
	std::vector<uint256_t> gc = compute_gc_dp_256(w, k);

	std::unordered_map<uint64_t, uint256_t> seen_masks;



	using Node256_64_64 = std::pair<uint256_t, std::pair<uint64_t, uint64_t>>;

	// Min-heap on (gc, mask, useless_mask)
	std::priority_queue<
		Node256_64_64,
		std::vector<Node256_64_64>,
		MinHeapComparator256_64_64
	> minHeap;

	// ---- seed singletons ----
	seen_masks[1ULL << 0] = gc[0] + 1;
	minHeap.emplace(gc[0] + 1, std::make_pair(1ULL << 0,0));
	for (uint64_t prefix = 1; prefix < (1ULL << (k - 1)); ++prefix) {
		uint64_t m = 1ULL << prefix;
		uint256_t v = gc[prefix] + 2;
		seen_masks[m] = v;
		minHeap.emplace(v, std::make_pair(m, 0ULL ));
	}

	uint256_t iterations = 0;
	while (!minHeap.empty()) {
		// --- check time limit only every 10,000 iterations to save in checks---
		if (iterations % 10000 == 0) {
			auto elapsed = clock::now() - start_time;
			if (elapsed > time_limit) {
				std::cout << "[timed] Time limit reached after "
					<< std::chrono::duration_cast<std::chrono::minutes>(elapsed).count()
					<< " min, returning best known number of charged windows." << std::endl;
				// Return 0 as mask, top heap gc as current minimum
				auto node = minHeap.top();

				uint256_t gc_value = node.first;
				uint64_t mask = node.second.first;
				uint64_t useless_mask = node.second.second;
				return { 0ULL, gc_value };
			}
			// cleanup every 30,000 iterations
			if (iterations % clean_every_n == 0 && iterations > 0) {


				// clean the heap ////////////////////////////////////
				auto cleanup_start = std::chrono::high_resolution_clock::now();
				size_t old_size = minHeap.size();

				std::cout << "[timed] Cleaning heap... current size: " << old_size << std::endl;

				std::vector<Node256_64_64> valid;

				valid.reserve(old_size);

				// Pop everything and keep only valid entries
				while (!minHeap.empty()) {
					auto node = minHeap.top();
					minHeap.pop();

					uint256_t gc_value = node.first;
					uint64_t mask = node.second.first;
					uint64_t useless_mask = node.second.second;

					auto it = seen_masks.find(mask);
					if (it != seen_masks.end() && it->second == gc_value)
						valid.emplace_back(gc_value, std::make_pair(mask, useless_mask));
				}

				// Reinsert valid entries
				for (auto& p : valid)
					minHeap.emplace(p);


				auto cleanup_end = std::chrono::high_resolution_clock::now();
				double cleanup_time =
					std::chrono::duration<double>(cleanup_end - cleanup_start).count();

				std::cout << "[timed] Heap cleaned, kept " << valid.size()
					<< " / " << old_size << " entries. Took "
					<< std::fixed << std::setprecision(3)
					<< cleanup_time << " s." << std::endl;


				// clean the seen_masks map /////////////////////////////
				uint256_t lowest_gc_in_heap = minHeap.top().first;

				// remove all entries with gc < lowest_gc_in_heap from seen_masks
				cleanup_start = std::chrono::high_resolution_clock::now();
				size_t before_cleanup = seen_masks.size();
				for (auto it = seen_masks.begin(); it != seen_masks.end(); ) {
					if (it->second < lowest_gc_in_heap) {
						it = seen_masks.erase(it);
					}
					else {
						++it;
					}
				}
				cleanup_end = std::chrono::high_resolution_clock::now();
				cleanup_time =
					std::chrono::duration<double>(cleanup_end - cleanup_start).count();
				std::cout << "[timed] Seen masks cleaned, kept " << seen_masks.size()
					<< " / " << before_cleanup << " entries. Took "
					<< std::fixed << std::setprecision(3)
					<< cleanup_time << " s." << std::endl;


			}
			
		}


		auto node = minHeap.top();
		minHeap.pop();

		uint256_t gc_value = node.first;
		uint64_t mask = node.second.first;
		uint64_t useless_mask = node.second.second;
		auto it_seen = seen_masks.find(mask);
		if (it_seen == seen_masks.end() || it_seen->second != gc_value)
			continue;

		iterations++;
		if (iterations % 1000 == 0) {
			std::cout << "Iterations:\t" << iterations
				<< "\theap_size:\t" << minHeap.size()
				<< "\tcurrent minimum CW:\t" << gc_value << std::endl;


			//analyze_mask(mask, useless_mask, k, w); // todo: remove later
		}

		

		bool can_be_uhs = true;
		std::vector<std::pair<uint64_t, uint256_t>> candidates;
		uint64_t free_mask = 0ULL;
		bool has_free = false;

		for (uint64_t kmer = 0; kmer < n; ++kmer) {
			if ((mask >> kmer) & 1ULL) continue;

			uint256_t temp = gc_value;

			uint256_t special_case = 0;
			if (kmer == 0 || kmer == n - 1)
				special_case = 1;
			


			const uint64_t new_mask = mask | (1ULL << kmer);
			const uint256_t pre = prefix_charged_dp_256(mask, kmer, w, k);
			const uint256_t suf = suffix_charged_dp_256_v2(mask, kmer, w, k);
			bool is_useless = false;
			if(pre+suf == 0) {
				// check if kmer is useless
				if (is_kmer_useless_dp_256(w, k, mask, kmer)) {
					is_useless = true;
					useless_mask |= (1ULL << kmer);
				}
			}
			else {
				can_be_uhs = false;

				//check if its in the edge case NEW NEW NEW - then we can ignore
				// actually assume w>=2, then context has at least 3 k-mers, then we can still make the same argument without checks!
				if (special_case == 0 && is_kmer_appearing_once(mask, kmer, k, w) && k>4) {
					continue;
				}


			}


			const uint256_t delta = pre + suf-special_case;

			if (delta== 0) {
				free_mask |= (1ULL << kmer);
				has_free = true;
			}
			else {
				temp += delta;
				if (temp <= n_charged)
					candidates.emplace_back(new_mask, temp);
			}
		}
		if (has_free) {

			// compute merged mask
			uint64_t merged_mask = mask | free_mask;

			// update seen_masks
			auto it = seen_masks.find(merged_mask);
			if (it == seen_masks.end() || gc_value < it->second) {
				seen_masks[merged_mask] = gc_value;

				// INSERT INTO HEAP: (gc_value, (merged_mask, new_useless_mask))
				minHeap.emplace(gc_value, std::make_pair(merged_mask, useless_mask));
			}

		}

		else {
			for (auto& [new_mask, temp] : candidates) {
				auto it = seen_masks.find(new_mask);
				if (it == seen_masks.end() || temp < it->second) {
					seen_masks[new_mask] = temp;
					minHeap.emplace(temp, std::make_pair(new_mask, useless_mask));
				}
			}
		}

		if (can_be_uhs && is_uhs_dp_wrapper_256(mask, w, k)) {
			std::cout << "[timed] reached a UHS (mask):\n" << std::bitset<64>(mask)
				<< " cw=" << gc_value << std::endl;

			// --- remove any useless kmers if recorded ---
			uint64_t cleaned_mask = mask & ~useless_mask;
			std::cout << "[timed] removed useless kmers, new set (mask):\n" << std::bitset<64>(cleaned_mask) << std::endl;
			return { cleaned_mask, gc_value };
		}
	}
	return { 0ULL, uint256_t(0) };
}

std::pair<uint64_t, boost::multiprecision::uint256_t>
dij_exhaust_dp_256(uint32_t w, uint32_t k, boost::multiprecision::uint256_t n_charged)
{
	using boost::multiprecision::uint256_t;

	std::cout << "w:\t" << w << "\tk:\t" << k << std::endl;

	const uint64_t n = 1ULL << k;               // k ≤ 6
	const uint64_t kmer_mask = n - 1ULL;
	std::vector<uint256_t> gc = compute_gc_dp_256(w, k);

	// seen: mask -> gc(256)
	std::unordered_map<uint64_t, uint256_t> seen_masks;

	// min-heap on (gc(256), mask(64))
	std::priority_queue<
		std::pair<uint256_t, uint64_t>,
		std::vector<std::pair<uint256_t, uint64_t>>,
		MinHeapComparator256_64
	> minHeap;

	// ---- seed singletons (0..0 and prefixes starting with 0) ----
	seen_masks[1ULL << 0] = gc[0] + 1;                   // correction for 0..0
	minHeap.emplace(gc[0] + 1, 1ULL << 0);
	for (uint64_t prefix = 1; prefix < (1ULL << (k - 1)); ++prefix) {
		uint64_t m = 1ULL << prefix;
		uint256_t v = gc[prefix] + 2;                    // correction for general prefix
		seen_masks[m] = v;
		minHeap.emplace(v, m);
	}

	uint256_t iterations = 0;
	while (!minHeap.empty()) {
		auto [gc_value, mask] = minHeap.top();
		minHeap.pop();
		auto it_seen = seen_masks.find(mask);
		if (it_seen == seen_masks.end() || it_seen->second != gc_value)
			continue;

		iterations++;
		if (iterations % 1000 == 0) {
			std::cout << "Iterations:\t" << iterations
				<< "\theap_size:\t" << minHeap.size()
				<< "\tcurrent_gc:\t" << gc_value << std::endl;
		}

		bool can_be_uhs = true;
		std::vector<std::pair<uint64_t, uint256_t>> candidates;
		uint64_t free_mask = 0ULL;  // bits for kmers with pre+suf == 0
		bool has_free = false;

		for (uint64_t kmer = 0; kmer < n; ++kmer) {
			if ((mask >> kmer) & 1ULL) continue;

			uint256_t temp = gc_value;
			if (kmer == 0 || kmer == n - 1)
				temp -= 1;

			const uint64_t new_mask = mask | (1ULL << kmer);
			const uint256_t pre = prefix_charged_dp_256(mask, kmer, w, k);
			const uint256_t suf = suffix_charged_dp_256_v2(mask, kmer, w, k);
			//const uint256_t suf = suffix_charged_dp_256(new_mask, kmer, w, k);
			const uint256_t delta = pre + suf;

			if (delta == 0) {
				free_mask |= (1ULL << kmer);
				has_free = true;
			}
			else {
				if (delta > 0) can_be_uhs = false;
				temp += delta;
				if (temp <= n_charged)
					candidates.emplace_back(new_mask, temp);
			}
		}

		if (has_free) {
			uint64_t merged_mask = mask | free_mask;
			// Keep same gc_value; only push if it's better
			auto it = seen_masks.find(merged_mask);
			if (it == seen_masks.end() || gc_value < it->second) {
				seen_masks[merged_mask] = gc_value;
				minHeap.emplace(gc_value, merged_mask);
			}
		}
		else {
			for (auto& [new_mask, temp] : candidates) {
				auto it = seen_masks.find(new_mask);
				if (it == seen_masks.end() || temp < it->second) {
					seen_masks[new_mask] = temp;
					minHeap.emplace(temp, new_mask);
				}
			}
		}

		if (can_be_uhs && is_uhs_dp_wrapper_256(mask, w, k)) {
			std::cout << "UHS found mask=" << mask << " gc=" << gc_value << std::endl;
			return { mask, gc_value };
		}
	}
	return { 0ULL, uint256_t(0) };
}


//std::pair<uint64_t, boost::multiprecision::uint256_t>
//dij_exhaust_parallel_shared_heap_dp_256(uint32_t w, uint32_t k,
//	boost::multiprecision::uint256_t n_charged, int N)
//{
//	using boost::multiprecision::uint256_t;
//
//	if (N < 1) N = 1;
//	const uint64_t n = 1ULL << k;
//	const uint64_t kmer_mask = n - 1ULL;
//
//	std::vector<uint256_t> gc = compute_gc_dp_256(w, k);
//
//	// sharded seen: mask(64) -> gc(256)
//	std::array<SeenShard64x256, N_SHARDS> seen_shards;
//	for (auto& s : seen_shards) { std::lock_guard<std::mutex> g(s.lock); s.map.reserve(200'000); }
//
//	// shared min-heap (gc(256), mask(64))
//	std::vector<std::pair<uint256_t, uint64_t>> heap;
//	heap.reserve(1'000'000);
//	std::mutex heap_mtx;
//	std::condition_variable heap_cv;
//
//	std::atomic<bool> stop{ false };
//	std::atomic<uint256_t> best_gc{ std::numeric_limits<uint64_t>::max() }; // large init ok
//	std::atomic<uint64_t> best_mask{ 0 };
//	std::atomic<bool> found{ false };
//	std::atomic<int>  active{ 0 };
//
//	// ---- Seed singletons ----
//	auto seed_push = [&](uint64_t m, const uint256_t& v) {
//		seen_set_64x256(seen_shards, m, v);
//		heap.emplace_back(v, m);
//		};
//
//	seed_push(1ULL << 0, gc[0] + 1);
//	for (uint64_t prefix = 1; prefix < (1ULL << (k - 1)); ++prefix) {
//		seed_push(1ULL << prefix, gc[prefix] + 2);
//	}
//	std::make_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//
//	// ---- Worker threads ----
//	auto worker = [&](int /*tid*/) {
//		std::vector<std::pair<uint256_t, uint64_t>> local;
//		local.reserve(MAX_LOCAL_PUSH_BATCH);
//
//		while (!stop) {
//			std::vector<std::pair<uint256_t, uint64_t>> batch;
//			{
//				std::unique_lock<std::mutex> lk(heap_mtx);
//				heap_cv.wait(lk, [&] { return stop || !heap.empty(); });
//				if (stop) return;
//
//				for (int i = 0; i < static_cast<int>(WORKLOAD_PER_WORKER) && !heap.empty(); ++i) {
//					std::pop_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//					batch.push_back(heap.back());
//					heap.pop_back();
//				}
//			}
//
//			active++;
//
//			for (auto& [gc_val, mask] : batch) {
//				uint256_t cur_seen = 0;
//				if (seen_lookup_64x256(seen_shards, mask, cur_seen) && cur_seen != gc_val)
//					continue;
//
//				bool can_be_uhs = true;
//
//				for (uint64_t kmer = 0; kmer < n; ++kmer) {
//					if ((mask >> kmer) & 1ULL) continue;
//
//					uint256_t tmp = gc_val;
//					if (kmer == 0 || kmer == (n - 1)) tmp -= 1;
//
//					const uint64_t new_mask = mask | (1ULL << kmer);
//					const uint256_t pre = prefix_charged_dp_256(mask, kmer, w, k);
//					const uint256_t suf = suffix_charged_dp_256(new_mask, kmer, w, k);
//					tmp += pre + suf;
//
//					if (tmp > n_charged) { can_be_uhs = false; continue; }
//					if (pre > 0 || suf > 0) can_be_uhs = false;
//
//					if (seen_update_if_better_64x256(seen_shards, new_mask, tmp)) {
//						local.emplace_back(tmp, new_mask);
//						if (local.size() >= MAX_LOCAL_PUSH_BATCH) {
//							std::lock_guard<std::mutex> lk(heap_mtx);
//							for (auto& p : local) {
//								heap.push_back(p);
//								std::push_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//							}
//							local.clear();
//							heap_cv.notify_all();
//						}
//					}
//				}
//
//				if (can_be_uhs && is_uhs_dp_wrapper(mask, w, k)) {
//					best_gc.store(gc_val);
//					best_mask.store(mask);
//					found.store(true);
//					stop.store(true);
//					heap_cv.notify_all();
//					break;
//				}
//			}
//
//			if (!local.empty()) {
//				std::lock_guard<std::mutex> lk(heap_mtx);
//				for (auto& p : local) {
//					heap.push_back(p);
//					std::push_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//				}
//				local.clear();
//				heap_cv.notify_all();
//			}
//
//			active--;
//		}
//		};
//
//	std::vector<std::thread> pool;
//	for (int i = 0; i < N; ++i) pool.emplace_back(worker, i);
//
//	// ---- Periodic cleanup ----
//	while (!stop) {
//		std::this_thread::sleep_for(std::chrono::milliseconds(CLEANUP_PERIOD_MS));
//		{
//			std::lock_guard<std::mutex> lk(heap_mtx);
//			rebuild_minheap_dedup_256_64(heap, seen_shards);
//			std::cout << "[256] cleanup, heap=" << heap.size() << std::endl;
//			if (heap.empty() && active.load() == 0) stop.store(true);
//		}
//		heap_cv.notify_all();
//	}
//
//	for (auto& t : pool) t.join();
//
//	if (found.load()) {
//		std::cout << "[256] Found UHS: mask=" << best_mask.load()
//			<< " gc=" << best_gc.load() << std::endl;
//		return { best_mask.load(), best_gc.load() };
//	}
//	return { 0ULL, uint256_t(0) };
//}
//
//
//
//std::pair<uint64_t, boost::multiprecision::uint256_t>
//dij_exhaust_parallel_shared_heap_dp_restrict_masks_256(
//	uint32_t w, uint32_t k,
//	boost::multiprecision::uint256_t n_charged,
//	int N)
//{
//	using boost::multiprecision::uint256_t;
//	if (N < 1) N = 1;
//	const uint64_t n = 1ULL << k;
//	const uint64_t kmer_mask = n - 1ULL;
//
//	std::vector<uint256_t> gc = compute_gc_dp_256(w, k);
//	std::array<SeenShard64x256, N_SHARDS> seen_shards;
//
//	std::vector<std::pair<uint256_t, uint64_t>> heap;
//	heap.reserve(1'000'000);
//	std::mutex heap_mutex;
//	std::condition_variable heap_cv;
//
//	std::atomic<bool> stop{ false };
//	std::atomic<uint256_t> best_gc{ std::numeric_limits<uint64_t>::max() };
//	std::atomic<uint64_t> best_mask{ 0 };
//	std::atomic<bool> uhs_found{ false };
//	std::atomic<int> active{ 0 };
//
//	// ---- Seed ----
//	auto seed_push = [&](uint64_t m, const uint256_t& v) {
//		seen_set_64x256(seen_shards, m, v);
//		heap.emplace_back(v, m);
//		};
//	seed_push(1ULL << 0, gc[0] + 1);
//	for (uint64_t prefix = 1; prefix < (1ULL << (k - 1)); ++prefix)
//		seed_push(1ULL << prefix, gc[prefix] + 2);
//	std::make_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//
//	// ---- Worker ----
//	auto worker = [&](int)
//		{
//			std::vector<std::pair<uint256_t, uint64_t>> local;
//			local.reserve(MAX_LOCAL_PUSH_BATCH);
//			std::vector<std::pair<uint256_t, uint64_t>> batch;
//			batch.reserve(WORKLOAD_PER_WORKER);
//
//			while (!stop) {
//				{
//					std::unique_lock<std::mutex> lk(heap_mutex);
//					heap_cv.wait(lk, [&] { return stop || !heap.empty(); });
//					if (stop) return;
//					for (int i = 0; i < WORKLOAD_PER_WORKER && !heap.empty(); ++i) {
//						std::pop_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//						batch.push_back(heap.back());
//						heap.pop_back();
//					}
//				}
//
//				active++;
//				for (auto& [gc_value, mask] : batch) {
//					uint256_t cur_seen = 0;
//					if (seen_lookup_64x256(seen_shards, mask, cur_seen) && cur_seen != gc_value)
//						continue;
//
//					bool can_be_uhs = true;
//					for (uint64_t kmer = 0; kmer < n; ++kmer) {
//						if ((mask >> kmer) & 1ULL) continue;
//						uint256_t tmp = gc_value;
//						if (kmer == 0 || kmer == (n - 1)) tmp -= 1;
//
//						const uint64_t new_mask = mask | (1ULL << kmer);
//						const uint256_t prefix_to_add = prefix_charged_dp_256(mask, kmer, w, k);
//						const uint256_t suffix_to_add = suffix_charged_dp_256(new_mask, kmer, w, k);
//
//						if (prefix_to_add == 0 && suffix_to_add == 0 &&
//							is_kmer_useless_dp(w, k, mask, kmer))
//							continue;
//
//						tmp += prefix_to_add + suffix_to_add;
//						if (tmp > n_charged) { can_be_uhs = false; continue; }
//						if (prefix_to_add > 0 || suffix_to_add > 0) can_be_uhs = false;
//
//						if (seen_update_if_better_64x256(seen_shards, new_mask, tmp)) {
//							local.emplace_back(tmp, new_mask);
//							if (local.size() >= MAX_LOCAL_PUSH_BATCH) {
//								std::lock_guard<std::mutex> lk(heap_mutex);
//								for (auto& p : local) {
//									heap.push_back(p);
//									std::push_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//								}
//								local.clear();
//								heap_cv.notify_all();
//							}
//						}
//					}
//
//					if (can_be_uhs && is_uhs_dp_wrapper(mask, w, k)) {
//						if (gc_value < best_gc.load()) {
//							best_gc.store(gc_value);
//							best_mask.store(mask);
//						}
//						uhs_found.store(true);
//						stop.store(true);
//						heap_cv.notify_all();
//					}
//				}
//
//				if (!local.empty()) {
//					std::lock_guard<std::mutex> lk(heap_mutex);
//					for (auto& p : local) {
//						heap.push_back(p);
//						std::push_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//					}
//					local.clear();
//					heap_cv.notify_all();
//				}
//
//				batch.clear();
//				active--;
//			}
//		};
//
//	std::vector<std::thread> pool;
//	for (int i = 0; i < N; ++i) pool.emplace_back(worker, i);
//
//	while (!stop) {
//		std::this_thread::sleep_for(std::chrono::milliseconds(CLEANUP_PERIOD_MS));
//		{
//			std::lock_guard<std::mutex> lk(heap_mutex);
//			rebuild_minheap_dedup_256_64(heap, seen_shards);
//			std::cout << "[256 restrict] heap=" << heap.size() << std::endl;
//			if (heap.empty() && active.load() == 0) stop.store(true);
//		}
//		heap_cv.notify_all();
//	}
//
//	for (auto& t : pool) t.join();
//
//	if (uhs_found) {
//		process_remaining_heap_single_dp_restrict_masks_256(
//			heap, seen_shards, w, k, n_charged, best_gc, best_mask);
//	}
//
//	if (best_mask.load() != 0) {
//		std::cout << "[256 restrict] Final mask=" << best_mask.load()
//			<< " gc=" << best_gc.load() << std::endl;
//		return { best_mask.load(), best_gc.load() };
//	}
//	return { 0ULL, uint256_t(0) };
//}

//static void process_remaining_heap_single_dp_256(
//std::vector<std::pair<boost::multiprecision::uint256_t, uint64_t>>& heap,
//std::array<SeenShard64x256, N_SHARDS>& seen_shards,
//uint32_t w, uint32_t k,
//boost::multiprecision::uint256_t n_charged,
//std::atomic<boost::multiprecision::uint256_t>& best_gc,
//std::atomic<uint64_t>& best_mask)
//{
//using boost::multiprecision::uint256_t;
//const uint64_t n = 1ULL << k;
//std::make_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//
//while (!heap.empty()) {
//	std::pop_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//	auto [gc_value, mask] = heap.back();
//	heap.pop_back();
//	if (gc_value >= best_gc.load()) break;
//
//	bool can_be_uhs = true;
//	for (uint64_t kmer = 0; kmer < n; ++kmer) {
//		if ((mask >> kmer) & 1ULL) continue;
//		uint256_t tmp = gc_value;
//		if (kmer == 0 || kmer == n - 1) tmp -= 1;
//
//		const uint64_t new_mask = mask | (1ULL << kmer);
//		const uint256_t pre = prefix_charged_dp_256(mask, kmer, w, k);
//		const uint256_t suf = suffix_charged_dp_256(new_mask, kmer, w, k);
//		tmp += pre + suf;
//		if (tmp > n_charged) { can_be_uhs = false; continue; }
//		if (pre > 0 || suf > 0) can_be_uhs = false;
//	}
//
//	if (can_be_uhs && is_uhs_dp_wrapper(mask, w, k)) {
//		if (gc_value < best_gc.load()) {
//			best_gc.store(gc_value);
//			best_mask.store(mask);
//			std::cout << "[Post-check 256] Better mask " << mask << " gc=" << gc_value << std::endl;
//		}
//	}
//}
//}

//static void process_remaining_heap_single_dp_restrict_masks_256(
//std::vector<std::pair<boost::multiprecision::uint256_t, uint64_t>>& heap,
//std::array<SeenShard64x256, N_SHARDS>& seen_shards,
//uint32_t w, uint32_t k,
//boost::multiprecision::uint256_t n_charged,
//std::atomic<boost::multiprecision::uint256_t>& best_gc,
//std::atomic<uint64_t>& best_mask)
//{
//using boost::multiprecision::uint256_t;
//const uint64_t n = 1ULL << k;
//std::make_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//
//while (!heap.empty()) {
//	std::pop_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//	auto [gc_value, mask] = heap.back();
//	heap.pop_back();
//	if (gc_value >= best_gc.load()) break;
//
//	uint256_t cur_seen = 0;
//	if (seen_lookup_64x256(seen_shards, mask, cur_seen) && cur_seen != gc_value)
//		continue;
//
//	bool can_be_uhs = true;
//	for (uint64_t kmer = 0; kmer < n; ++kmer) {
//		if ((mask >> kmer) & 1ULL) continue;
//
//		uint256_t tmp = gc_value;
//		if (kmer == 0 || kmer == n - 1) tmp -= 1;
//
//		const uint64_t new_mask = mask | (1ULL << kmer);
//		const uint256_t pre = prefix_charged_dp_256(mask, kmer, w, k);
//		const uint256_t suf = suffix_charged_dp_256(new_mask, kmer, w, k);
//
//		if (pre == 0 && suf == 0 && is_kmer_useless_dp(w, k, mask, kmer))
//			continue;
//
//		tmp += pre + suf;
//		if (tmp > n_charged) { can_be_uhs = false; continue; }
//		if (pre > 0 || suf > 0) can_be_uhs = false;
//
//		if (seen_update_if_better_64x256(seen_shards, new_mask, tmp) && tmp < best_gc.load()) {
//			heap.emplace_back(tmp, new_mask);
//			std::push_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//		}
//	}
//
//	if (can_be_uhs && is_uhs_dp_wrapper(mask, w, k)) {
//		auto old_best = best_gc.load();
//		if (gc_value < old_best) {
//			best_gc.store(gc_value);
//			best_mask.store(mask);
//			std::cout << "[Post-check 256 restrict] Better gc=" << gc_value << " mask=" << mask << std::endl;
//
//			heap.erase(std::remove_if(heap.begin(), heap.end(),
//				[&](const auto& p) { return p.first >= gc_value; }),
//				heap.end());
//			std::make_heap(heap.begin(), heap.end(), MinHeapComparator256_64{});
//		}
//	}
//}
//}
//

std::vector<uint64_t> get_best_order_dp_256(
	uint64_t original_mask,
	boost::multiprecision::uint256_t best_gc,
	uint32_t w, uint32_t k)
{
	using boost::multiprecision::uint256_t;
	const uint64_t n = 1ULL << k;
	std::vector<uint256_t> gc = compute_gc_dp_256(w, k);

	std::unordered_map<uint64_t, uint256_t> old_mask_to_gc;
	std::unordered_map<uint64_t, std::vector<uint64_t>> mask_to_order;

	// ---- initialize ----
	if (original_mask & 1ULL) {
		old_mask_to_gc[1ULL] = gc[0] + 1;
		mask_to_order[1ULL] = { 0 };
	}
	for (uint64_t prefix = 1; prefix < (1ULL << (k - 1)); ++prefix) {
		if ((original_mask >> prefix) & 1ULL) {
			uint64_t m = 1ULL << prefix;
			old_mask_to_gc[m] = gc[prefix] + 2;
			mask_to_order[m] = { prefix };
		}
	}

	// ---- iterative reconstruction ----
	for (uint64_t rank = 0; rank <= n; ++rank) {
		std::unordered_map<uint64_t, uint256_t> mask_to_gc;

		for (auto& [mask, value] : old_mask_to_gc) {
			if (is_uhs_dp_wrapper_256(mask, w, k)) {
				std::cout << "UHS found, number of charged windows:\n" << value << std::endl;
				return mask_to_order[mask];
			}

			uint64_t free_mask = 0ULL;
			bool has_free = false;
			std::vector<std::pair<uint64_t, uint256_t>> candidates;

			for (uint64_t kmer = 0; kmer < n; ++kmer) {
				// skip kmers not in the original
				if (!((original_mask >> kmer) & 1ULL))
					continue;

				if ((mask >> kmer) & 1ULL)
					continue;

				uint256_t temp_value = value;

				uint256_t special_case = 0;
				if (kmer == 0 || kmer == n - 1)
					special_case = 1;

				const uint64_t new_mask = mask | (1ULL << kmer);
				const uint256_t pre = prefix_charged_dp_256(mask, kmer, w, k);
				//const uint256_t suf = suffix_charged_dp_256(new_mask, kmer, w, k);
				const uint256_t suf = suffix_charged_dp_256_v2(mask, kmer, w, k);

				const uint256_t delta = pre + suf- special_case;

				if (delta == 0) {
					if (!is_kmer_useless_dp_256(w, k, mask, kmer)) {
						free_mask |= (1ULL << kmer);
						has_free = true;
					}
					else {
						continue;
					}
				}
				else {
					temp_value += delta;
					if (temp_value <= best_gc)
						candidates.emplace_back(new_mask, temp_value);
				}
			}

			if (has_free) {
				uint64_t merged_mask = mask | free_mask;
				if (mask_to_gc.find(merged_mask) == mask_to_gc.end() ||
					value < mask_to_gc[merged_mask]) {
					mask_to_gc[merged_mask] = value;
					mask_to_order[merged_mask] = mask_to_order[mask];
					for (uint64_t kmer = 0; kmer < n; ++kmer) {
						if ((free_mask >> kmer) & 1ULL)
							mask_to_order[merged_mask].push_back(kmer);
					}
				}
			}
			else {
				for (auto& [new_mask, temp_value] : candidates) {
					if (mask_to_gc.find(new_mask) == mask_to_gc.end() ||
						temp_value < mask_to_gc[new_mask]) {
						mask_to_gc[new_mask] = temp_value;
						mask_to_order[new_mask] = mask_to_order[mask];
						// append kmer that led to this mask
						for (uint64_t kmer = 0; kmer < n; ++kmer)
							if (((new_mask ^ mask) >> kmer) & 1ULL)
								mask_to_order[new_mask].push_back(kmer);
					}
				}
			}
		}

		old_mask_to_gc = std::move(mask_to_gc);
		if (old_mask_to_gc.empty()) {
			std::cout << "[256] Cannot reconstruct order" << std::endl;
			return {};
		}
	}
	return {};
}



std::pair<std::vector<uint64_t>,uint256_t> get_best_order_and_gc_dp_256(
	uint64_t original_mask,
	boost::multiprecision::uint256_t best_gc,
	uint32_t w, uint32_t k, bool verbose=true)
{
	using boost::multiprecision::uint256_t;
	const uint64_t n = 1ULL << k;
	std::vector<uint256_t> gc = compute_gc_dp_256(w, k);

	std::unordered_map<uint64_t, uint256_t> old_mask_to_gc;
	std::unordered_map<uint64_t, std::vector<uint64_t>> mask_to_order;

	// ---- initialize ----
	if (original_mask & 1ULL) {
		old_mask_to_gc[1ULL] = gc[0] + 1;
		mask_to_order[1ULL] = { 0 };
	}
	for (uint64_t prefix = 1; prefix < (1ULL << (k - 1)); ++prefix) {
		if ((original_mask >> prefix) & 1ULL) {
			uint64_t m = 1ULL << prefix;
			old_mask_to_gc[m] = gc[prefix] + 2;
			mask_to_order[m] = { prefix };
		}
	}

	// ---- iterative reconstruction ----
	for (uint64_t rank = 0; rank <= n; ++rank) {
		std::unordered_map<uint64_t, uint256_t> mask_to_gc;

		for (auto& [mask, value] : old_mask_to_gc) {
			if (is_uhs_dp_wrapper_256(mask, w, k)) {
				if (verbose) {
					std::cout << "UHS found, gc=" << value << std::endl;
				}
                return std::make_pair(mask_to_order[mask], value);
			}

			uint64_t free_mask = 0ULL;
			bool has_free = false;
			std::vector<std::pair<uint64_t, uint256_t>> candidates;

			for (uint64_t kmer = 0; kmer < n; ++kmer) {
				// skip kmers not in the original
				if (!((original_mask >> kmer) & 1ULL))
					continue;

				if ((mask >> kmer) & 1ULL)
					continue;

				uint256_t temp_value = value;
				if (kmer == 0 || kmer == n - 1)
					temp_value -= 1;

				const uint64_t new_mask = mask | (1ULL << kmer);
				const uint256_t pre = prefix_charged_dp_256(mask, kmer, w, k);
				//const uint256_t suf = suffix_charged_dp_256(new_mask, kmer, w, k);
				const uint256_t suf = suffix_charged_dp_256_v2(mask, kmer, w, k);

				const uint256_t delta = pre + suf;

				if (delta == 0) {
					free_mask |= (1ULL << kmer);
					has_free = true;
				}
				else {
					temp_value += delta;
					if (temp_value <= best_gc)
						candidates.emplace_back(new_mask, temp_value);
				}
			}

			if (has_free) {
				uint64_t merged_mask = mask | free_mask;
				if (mask_to_gc.find(merged_mask) == mask_to_gc.end() ||
					value < mask_to_gc[merged_mask]) {
					mask_to_gc[merged_mask] = value;
					mask_to_order[merged_mask] = mask_to_order[mask];
					for (uint64_t kmer = 0; kmer < n; ++kmer) {
						if ((free_mask >> kmer) & 1ULL)
							mask_to_order[merged_mask].push_back(kmer);
					}
				}
			}
			else {
				for (auto& [new_mask, temp_value] : candidates) {
					if (mask_to_gc.find(new_mask) == mask_to_gc.end() ||
						temp_value < mask_to_gc[new_mask]) {
						mask_to_gc[new_mask] = temp_value;
						mask_to_order[new_mask] = mask_to_order[mask];
						// append kmer that led to this mask
						for (uint64_t kmer = 0; kmer < n; ++kmer)
							if (((new_mask ^ mask) >> kmer) & 1ULL)
								mask_to_order[new_mask].push_back(kmer);
					}
				}
			}
		}

		old_mask_to_gc = std::move(mask_to_gc);
		if (old_mask_to_gc.empty()) {
			std::cout << "[256] Cannot reconstruct order" << std::endl;
			return {};
		}
	}
	return {};
}
