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
#include <chrono>
#include <thread>
#include <future>
#include <omp.h>
#include <iomanip>  // for std::setprecision
#include <fstream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <queue>
#include <atomic>
#include <set>

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







using boost::multiprecision::uint256_t;



void paper_k_six_lower_bound_parallel()
{
	using boost::multiprecision::uint256_t;

	const uint32_t k = 6;
	const double time_limit_minutes = 60.0*4;
	//const double time_limit_minutes = 10;

	const int MAX_THREADS = 16;

	// --- Define w values to test ---
	// vector is FIFO
	std::vector<uint32_t> w_values = {};
	//for (uint32_t w = 101; w <= 200; ++w) // should be 200
	//	w_values.insert(w);
	// Single values
	//w_values.insert(138);


	// Ranges
	//for (uint32_t w = 148; w <= 157; ++w)
	//	w_values.insert(w);

	//for (uint32_t w = 20; w <= 118; ++w)
	//	w_values.push_back(w);
	// 


	for (uint32_t w = 36; w >= 2; --w)
		w_values.push_back(w);

	// 
	// TODO NEXT
	//for (uint32_t w = 125; w <= 129; ++w)
	//	w_values.insert(w);




	// --- CSV output setup ---
	std::ofstream csv("paper_k6_lower_bounds_parallel.csv");
	csv << "w,k,best_gc,best_mask,new_order,best_gc_lower_bound,time_seconds\n";
	std::mutex csv_mutex;

	// --- Work queue ---
	std::queue<uint32_t> work_queue;
	for (auto w : w_values) work_queue.push(w);
	std::mutex queue_mutex;

	std::atomic<int> active_threads(0);

	auto worker = [&]() {
		while (true) {
			uint32_t w;
			{   // lock scope
				std::lock_guard<std::mutex> lock(queue_mutex);
				if (work_queue.empty()) break;
				w = work_queue.front();
				work_queue.pop();
			}

			active_threads++;
			std::cout << "\n==============================\n";
			std::cout << "Thread " << std::this_thread::get_id()
				<< " running w=" << w << ", k=" << k << std::endl;

			auto start = std::chrono::high_resolution_clock::now();

			auto result = dij_exhaust_dp_256_time_limit(
				w, k, std::numeric_limits<uint256_t>::max(), time_limit_minutes, 30000);

			uint64_t best_mask = result.first;
			uint256_t best_gc = result.second;

			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = end - start;
			double time_sec = elapsed.count();

			std::ostringstream line;

			if (best_mask != 0) {
				std::cout << "[Finished] w=" << w
					<< " best_gc=" << best_gc
					<< " mask=" << std::bitset<64>(best_mask) << std::endl;

				std::vector<uint64_t> order = get_best_order_dp_256(best_mask, best_gc, w, k);
				std::vector<uint64_t> new_order = convert_to_order(order, k);

				std::ostringstream new_order_str;
				new_order_str << "\"";
				for (size_t i = 0; i < new_order.size(); ++i) {
					new_order_str << new_order[i];
					if (i + 1 < new_order.size()) new_order_str << " ";
				}
				new_order_str << "\"";

				line << w << "," << k << ","
					<< best_gc << ","
					<< best_mask << ","
					<< new_order_str.str() << ","
					<< ","  // empty lower bound column
					<< std::fixed << std::setprecision(3) << time_sec << "\n";
			}
			else {
				std::cout << "[Timeout] w=" << w
					<< " lower_bound_gc=" << best_gc
					<< " (limit " << time_limit_minutes << " min)" << std::endl;

				line << w << "," << k << ","
					<< ","  // best_gc
					<< ","  // best_mask
					<< ","  // new_order
					<< best_gc << ","
					<< std::fixed << std::setprecision(3) << time_sec << "\n";
			}

			{   // safely write to CSV
				std::lock_guard<std::mutex> lock(csv_mutex);
				csv << line.str();
				csv.flush();
			}

			active_threads--;
		}
		};

	// --- Launch up to MAX_THREADS threads ---
	std::vector<std::thread> threads;
	for (int i = 0; i < MAX_THREADS; ++i)
		threads.emplace_back(worker);

	// --- Wait for all threads ---
	for (auto& t : threads) t.join();

	csv.close();
	std::cout << "\n✅ Results saved to paper_k6_lower_bounds_parallel.csv\n";
}

void paper_k_small_parallel()
{
	using boost::multiprecision::uint256_t;

	const double time_limit_minutes = 60.0;
	const int MAX_THREADS = 1;

	// --- Build work list for (k,w) pairs ---
	std::queue<std::pair<uint32_t, uint32_t>> work_queue;


	for (uint32_t k = 5; k <= 5; ++k) {
		uint32_t w_start = 2;
		uint32_t w_end = 3 * (1U << k); // 3 * 2^k
		for (uint32_t w = w_start; w <= w_end; ++w)
			work_queue.emplace(k, w);
	}

	//for (uint32_t k = 2; k <= 5; ++k) {
	//	uint32_t w_start = 6;
	//	uint32_t w_end = 3 * (1U << k); // 3 * 2^k
	//	for (uint32_t w = w_start; w <= w_end; ++w)
	//		work_queue.emplace(k, w);
	//}


	// --- CSV output setup ---
	std::ofstream csv("paper_k2_to_k5_lower_bounds_parallel.csv");
	csv << "w,k,best_gc,best_mask,new_order,best_gc_lower_bound,time_seconds\n";
	std::mutex csv_mutex;

	std::mutex queue_mutex;
	std::atomic<int> active_threads(0);

	auto worker = [&]() {
		while (true) {
			uint32_t k, w;
			{   // lock scope for work queue
				std::lock_guard<std::mutex> lock(queue_mutex);
				if (work_queue.empty()) break;
				k = work_queue.front().first;
				w = work_queue.front().second;
				work_queue.pop();
			}

			active_threads++;
			std::cout << "\n==============================\n";
			std::cout << "Thread " << std::this_thread::get_id()
				<< " running w=" << w << ", k=" << k << std::endl;

			auto start = std::chrono::high_resolution_clock::now();

			auto result = dij_exhaust_dp_256_time_limit(
				w, k, std::numeric_limits<uint256_t>::max(), time_limit_minutes, 1000000);

			uint64_t best_mask = result.first;
			uint256_t best_gc = result.second;

			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = end - start;
			double time_sec = elapsed.count();

			std::ostringstream line;

			if (best_mask != 0) {
				std::cout << "[Finished] w=" << w << " k=" << k
					<< " best_gc=" << best_gc
					<< " mask=" << std::bitset<64>(best_mask) << std::endl;

				std::vector<uint64_t> order = get_best_order_dp_256(best_mask, best_gc, w, k);
				std::vector<uint64_t> new_order = convert_to_order(order, k);

				std::ostringstream new_order_str;
				new_order_str << "\"";
				for (size_t i = 0; i < new_order.size(); ++i) {
					new_order_str << new_order[i];
					if (i + 1 < new_order.size()) new_order_str << " ";
				}
				new_order_str << "\"";

				line << w << "," << k << ","
					<< best_gc << ","
					<< best_mask << ","
					<< new_order_str.str() << ","
					<< ","  // empty lower bound column
					<< std::fixed << std::setprecision(3) << time_sec << "\n";
			}
			else {
				std::cout << "[Timeout] w=" << w << " k=" << k
					<< " lower_bound_gc=" << best_gc
					<< " (limit " << time_limit_minutes << " min)" << std::endl;

				line << w << "," << k << ","
					<< ","  // best_gc
					<< ","  // best_mask
					<< ","  // new_order
					<< best_gc << ","
					<< std::fixed << std::setprecision(3) << time_sec << "\n";
			}

			{   // safely write to CSV
				std::lock_guard<std::mutex> lock(csv_mutex);
				csv << line.str();
				csv.flush();
			}

			active_threads--;
		}
		};

	// --- Launch up to MAX_THREADS threads ---
	std::vector<std::thread> threads;
	for (int i = 0; i < MAX_THREADS; ++i)
		threads.emplace_back(worker);

	// --- Wait for all threads ---
	for (auto& t : threads)
		t.join();

	csv.close();
	std::cout << "\n✅ Results saved to paper_k2_to_k5_lower_bounds_parallel.csv\n";
}



void paper_k_six_lower_bound_routine() {
	using boost::multiprecision::uint256_t;

	std::ofstream csv("paper_k6_lower_bounds.csv");
	csv << "w,k,best_gc,best_mask,new_order,best_gc_lower_bound,time_seconds\n";

	uint32_t k = 6;
	double time_limit_minutes = 20.0; 

	for (uint32_t w = 65; w <= 205; w+=10) {
		std::cout << "\n==============================" << std::endl;
		std::cout << "Running w=" << w << ", k=" << k << std::endl;
		auto start = std::chrono::high_resolution_clock::now();

		auto result = dij_exhaust_dp_256_time_limit(w, k, std::numeric_limits<uint256_t>::max(), time_limit_minutes, 30000);
		uint64_t best_mask = result.first;
		uint256_t best_gc = result.second;

		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end - start;
		double time_sec = elapsed.count();

		std::ostringstream new_order_str;

		// Check if completed (mask != 0)
		if (best_mask != 0) {
			std::cout << "[Finished] w=" << w << " k=" << k
				<< " best_gc=" << best_gc
				<< " mask=" << std::bitset<64>(best_mask) << std::endl;

			std::vector<uint64_t> order = get_best_order_dp_256(best_mask, best_gc, w, k);
			std::vector<uint64_t> new_order = convert_to_order(order, k);

			new_order_str << "\""; // CSV-safe string
			for (size_t i = 0; i < new_order.size(); ++i) {
				new_order_str << new_order[i];
				if (i + 1 < new_order.size()) new_order_str << " ";
			}
			new_order_str << "\"";

			csv << w << "," << k << ","
				<< best_gc << ","
				<< best_mask << ","
				<< new_order_str.str() << ","
				<< ","  // empty best_gc_lower_bound column
				<< std::fixed << std::setprecision(3) << time_sec << "\n";
		}
		else {
			std::cout << "[Timeout] w=" << w << " k=" << k
				<< " lower_bound_gc=" << best_gc
				<< " (limit " << time_limit_minutes << " min)" << std::endl;

			csv << w << "," << k << ","
				<< ","  // best_gc empty
				<< ","  // best_mask empty
				<< ","  // new_order empty
				<< best_gc << ","  // best_gc_lower_bound filled
				<< std::fixed << std::setprecision(3) << time_sec << "\n";
		}
	}

	csv.close();
	std::cout << "\nResults saved to paper_k6_lower_bounds.csv ✅" << std::endl;
}


void k_six_routine() {
	// start measuring time

    auto start = std::chrono::high_resolution_clock::now();

	// run for k=6 w=200

	

	uint32_t k = 6;
	uint32_t w = 192;

	std::cout << "running for k=6 w=192" << std::endl;

	std::pair<uint64_t, uint256_t> result = dij_exhaust_dp_256(w, k, std::numeric_limits<uint256_t>::max());
	uint64_t best_mask = result.first;
	uint256_t best_gc = result.second;

	std::cout << "w:\t" << w << "\tk:\t" << k << "\tbest_gc:\t" << best_gc << "\tbest_mask:\t" << std::bitset<64>(best_mask) << std::endl;

	std::vector<uint64_t> order = get_best_order_dp_256(best_mask, best_gc, w, k);


	std::vector<uint64_t> new_order = convert_to_order(order, k);
	std::cout << "best order:\t";
	for (auto& kmer : new_order) {
		std::cout << kmer << ", ";
	}
	std::cout << std::endl;

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Time taken by function: "
			  << elapsed.count() << " seconds" << std::endl;





}

void large_w_routine() {





	// run for 2<=k<=5
	for (uint32_t k = 2; k <= 5; k++) {


		for (uint32_t w = 5; w <= 56; w++) {

			std::cout << "\n\n\n";
			std::vector<uint64_t> order = {};


			if (k == 5) {
				std::pair<uint64_t, uint64_t> result = dij_exhaust_parallel_shared_heap_dp_restrict_masks(w, k, std::numeric_limits<uint64_t>::max());

				uint64_t best_mask = result.first;
				uint64_t best_gc = result.second;

				std::cout << "w:\t" << w << "\tk:\t" << k << "\tbest_gc:\t" << best_gc << "\tbest_mask:\t" << std::bitset<64>(best_mask) << std::endl;

				order = get_best_order_dp(best_mask, best_gc, w, k);
			}
			else {
				std::pair<uint64_t, uint64_t> result = {};
				if (w < 15) {
					result = dij_exhaust_dfs(w, k, std::numeric_limits<uint64_t>::max());

				}
				else {
					result = dij_exhaust_dp(w, k, std::numeric_limits<uint64_t>::max());
				}

				uint64_t best_mask = result.first;
				uint64_t best_gc = result.second;

				std::cout << "w:\t" << w << "\tk:\t" << k << "\tbest_gc:\t" << best_gc << "\tbest_mask:\t" << std::bitset<64>(best_mask) << std::endl;

				if (w > 15) {
					order = get_best_order_dp(best_mask, best_gc, w, k);

				}
				else {
					order = get_best_order(best_mask, best_gc, w, k);
				}

			}







			std::vector<uint64_t> new_order = convert_to_order(order, k);
			std::cout << "best order:\t";
			for (auto& kmer : new_order) {
				std::cout << kmer << ", ";
			}
			std::cout << std::endl;
		}



		for (uint32_t w = 57; w <= 100; w++) {

			std::cout << "\n\n\n";
			std::pair<uint64_t, uint256_t> result = dij_exhaust_dp_256(w, k, std::numeric_limits<uint256_t>::max());

			uint64_t best_mask = result.first;
			uint256_t best_gc = result.second;

			std::cout << "w:\t" << w << "\tk:\t" << k << "\tbest_gc:\t" << best_gc << "\tbest_mask:\t" << std::bitset<64>(best_mask) << std::endl;

			std::vector<uint64_t> order = get_best_order_dp_256(best_mask, best_gc, w, k);

			std::vector<uint64_t> new_order = convert_to_order(order, k);

			// print the order

			std::cout << "best order:\t";
			for (auto& kmer : new_order) {
				std::cout << kmer << ", ";
			}
			std::cout << std::endl;
		}
	}
}



void optimal_orders_routine(uint32_t w, uint32_t k, double time_limit_minutes = 60.0) {


	using boost::multiprecision::uint256_t;


	auto start = std::chrono::high_resolution_clock::now();


	uint64_t clean_every_n = 0;
	if (k == 6) {
		clean_every_n = 30000;
	}
	else {
		clean_every_n = 1000000;
	}

	auto result = dij_exhaust_dp_256_time_limit(
		w, k, std::numeric_limits<uint256_t>::max(), time_limit_minutes, clean_every_n);

	uint64_t best_mask = result.first;
	uint256_t best_gc = result.second;

	auto end = std::chrono::high_resolution_clock::now();


	std::ostringstream line;

	if (best_mask != 0) {
		std::cout << "Reconstructing order from set (mask) of k-mers..." << std::endl;

		std::vector<uint64_t> order = get_best_order_dp_256(best_mask, best_gc, w, k);
		std::vector<uint64_t> new_order = convert_to_order(order, k);

		std::chrono::duration<double> elapsed = end - start;
		double time_sec = elapsed.count();

		std::ostringstream new_order_str;
		new_order_str << "";
		for (size_t i = 0; i < new_order.size(); ++i) {
			new_order_str << new_order[i];
			if (i + 1 < new_order.size()) new_order_str << " ";
		}
		new_order_str << "";

		std::cout << "Order: \n" << new_order_str.str() << std::endl;
		std::cout << "Time taken: " << time_sec << " seconds" << std::endl;
	}
	else {
		std::cout << "[Timeout] w=" << w
			<< " Lower bound for number of charged windows:\n" << best_gc
			<< "\n(limit " << time_limit_minutes << " min)" << std::endl;
	}

}


void paper_optimal_uhs_routine()
{
	// Open CSV file once, write header
	std::ofstream csv("uhs_summary.csv");
	csv << "k,w,n_optimal_uhs,uhs_size_best,best_gc_uhs,best_gc_opt,uhs_size_opt\n";

	std::vector<std::pair<uint32_t, uint32_t>> duos;
	for (uint32_t k = 2; k <= 4; k++) {
		uint32_t w_start = 2;
		if (k == 4) {
			w_start = 4;
		}
		uint32_t w_end = (1ULL << k) * 3; // 3 * 2^k
		for (uint32_t w = w_start; w <= w_end; w++) {
			duos.emplace_back(k, w);
		}
	}

	for (auto [k, w] : duos) {
		std::cout << "\n\n k=" << k << " w=" << w << std::endl;
		auto [opt_uhs, opt_uhs_size] = get_opt_uhs_dp_256(w, k);

		std::list<uint256_t> all_gcs;
		uint256_t best_gc = std::numeric_limits<uint256_t>::max();
		uint64_t best_mask = 0;
		uint32_t counter = 0;
		uint256_t max_uint256 = std::numeric_limits<uint256_t>::max();

		for (uint64_t mask : opt_uhs) {
			auto [order, gc] = get_best_order_and_gc_dp_256(mask, max_uint256, w, k, false);
			if (gc < best_gc) {
				best_gc = gc;
				best_mask = mask;
			}
			all_gcs.push_back(gc);
			counter++;

			if (counter % 1000 == 0) {
				std::cout << "Checked " << counter
					<< " optimal UHS masks, best gc so far: " << best_gc << std::endl;
			}
		}

		uint64_t n_optimal_uhs = opt_uhs.size();
		std::cout << "Best UHS mask:\t" << std::bitset<64>(best_mask)
			<< "\twith gc:\t" << best_gc << std::endl;

		uint32_t uhs_size_best = std::popcount(best_mask);
		uint256_t best_gc_uhs = best_gc;


		// Compute optimal order by exhaustive DP
		auto [opt_mask, opt_gc] = dij_exhaust_dp_256_time_limit(w, k, max_uint256, 60,1000000);
		std::vector<uint64_t> opt_order = get_best_order_dp_256(opt_mask, opt_gc, w, k);
		std::cout << "Optimal gc for w=" << w << " k=" << k << " is " << opt_gc << std::endl;
		uint256_t best_gc_opt = opt_gc;

		uint32_t uhs_size_opt = get_uhs_size_dp_256(opt_order, w, k);
		std::cout << "Optimal order UHS size for w=" << w << " k=" << k
			<< " is " << uhs_size_opt << std::endl;

		// ---- ✅ Write single summary line to CSV ----
		csv << k << ","
			<< w << ","
			<< n_optimal_uhs << ","
			<< uhs_size_best << ","
			<< best_gc_uhs << ","
			<< best_gc_opt << ","
			<< uhs_size_opt << "\n";

		csv.flush(); // optional: ensure data written even if crash
	}

	csv.close();
}



void optimal_uhs_routine(uint32_t w, uint32_t k) {
	auto [opt_uhs, opt_uhs_size] = get_opt_uhs_dfs(w, k);



	// for each opt_uhs get opt order and keep the one with the lowest GC
	std::list<uint64_t> all_gcs = std::list<uint64_t>();
	uint64_t best_gc = std::numeric_limits<uint64_t>::max();
	uint64_t best_mask = 0;
	uint32_t counter = 0;
	for (uint64_t mask : opt_uhs) {
		auto [order, gc] = get_best_order_from_uhs(mask, opt_uhs_size, w, k);
		if (gc < best_gc) {
			best_gc = gc;
			best_mask = mask;
		}
		all_gcs.push_back(gc);
		counter++;

		if (counter % 100 == 0) {
			std::cout << "Checked\t" << counter << "\toptimal UHS masks, best gc so far:\t" << best_gc << std::endl;
		}
	}

	std::cout << "Best UHS mask:\t" << std::bitset<64>(best_mask) << "\twith gc:\t" << best_gc << std::endl;
	// save all the gcs to a list in a file
	std::ofstream outfile("all_uhs_gcs_w_" + std::to_string(w) + "_k_" + std::to_string(k) + ".txt");
	for (uint64_t gc : all_gcs) {
		outfile << gc << std::endl;
	}
}




void random_density_routine(uint32_t w, uint32_t k) {
	double density = get_random_density_256(w, k);
	std::cout << "Average density factor for w:\t" << w << "\tk:\t" << k << "\tis\t" << density << std::endl;

}


double run_with_timeout(uint32_t w, uint32_t k, std::chrono::minutes timeout) {
	// Run get_random_density_256 asynchronously
	auto future = std::async(std::launch::async, [w, k]() {
		return get_random_density_256(w, k);
		});

	// Wait for result or timeout
	if (future.wait_for(timeout) == std::future_status::ready) {
		return future.get();
	}
	else {
		std::cerr << "Timeout: k=" << k << " w=" << w << " exceeded "
			<< timeout.count() << " minutes.\n";
		return -1.0; // sentinel value
	}
}

void paper_random_density_routine() {
	// --- CSV for normal ranges ---
	std::ofstream csv("random_density_summary.csv");
	csv << "k,w,density\n";

	// Iterate k=2,3,4 and w=2..100
	for (uint32_t k = 2; k <= 4; ++k) {
		for (uint32_t w = 2; w <= 100; ++w) {
			double density = get_random_density_256(w, k);
			csv << k << "," << w << "," << density << "\n";
			std::cout << "k=" << k << " w=" << w
				<< " density=" << density << std::endl;
		}
	}

	csv.close();

	// --- Separate test for heavy (k=5, w=200) with 20 min limit ---
	std::ofstream csv_heavy("random_density_k5_w200.csv");
	csv_heavy << "k,w,density\n";

	uint32_t k = 5, w = 200;
	double density_heavy = run_with_timeout(w, k, std::chrono::minutes(20));
	csv_heavy << k << "," << w << "," << density_heavy << "\n";
	std::cout << "k=" << k << " w=" << w
		<< " density=" << density_heavy << " (20min limit)" << std::endl;

	csv_heavy.close();
}





int main(int  argc, char* argv[])
{
	// base things user runs
	//paper_optimal_uhs_routine();
	//random_density_routine(15,4);
	//optimal_orders_routine(192, 6, 10);
	//return 0;


	// things for testing
	//paper_k_small_parallel();
	//paper_k_six_lower_bound_parallel();
	//paper_optimal_uhs_routine();
	//paper_random_density_routine();
	//return 0;
	



	uint32_t w = 0, k = 0;
    std::string routine;
	double time_limit_minutes = 60.0;
	uint32_t n_cores = 1;


    // Simple manual argument parser
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-w" && i + 1 < argc) {
            w = static_cast<uint32_t>(std::stoi(argv[++i]));
        }
        else if (arg == "-k" && i + 1 < argc) {
            k = static_cast<uint32_t>(std::stoi(argv[++i]));
        }
        else if (arg == "-routine" && i + 1 < argc) {
            routine = argv[++i];
        }
		else if (arg == "-time_limit_minutes" && i + 1 < argc) {
			time_limit_minutes = std::stod(argv[++i]);
		}
        else {
            std::cerr << "Unknown or incomplete argument: " << arg << std::endl;
            return 1;
        }
    }

	///111814
	
	// FOR TESTING SET IT AS CONSTANTS
	//w = 20;
	//k = 4;
	//routine = "OptMini";
	
	//routine = "uhs_tests";
	//routine = "calc_average_dens";







	// validate if picked any routine
	if (routine.empty()) {
		std::cerr << "Error: -routine is a required parameter.\n"
			<< "Usage: ./my_program -routine <OM-heap|calc_average_dens> -w <w> -k <k> -time_limit_minutes <time_limit_in_minutes>\n";
	}





    // Select routine
    if (routine == "OM-heap") {

		// Validate required args
		if (w == 0 || k == 0) {
			std::cerr << "Error: -w and -k are required parameters.\n"
				<< "Usage: ./my_program -routine <OM-heap|calc_average_dens> -w <w> -k <k> -time_limit_minutes <time_limit_in_minutes>\n";
			return 1;
		}
		// Validate required args
		if (w < 2) {
			std::cerr << "Error: -w must be at least 2.\n";
			return 1;
		}
		if (k < 1) {
			std::cerr << "Error: -k must be at least 1.\n";
			return 1;
		}
		if (k > 6) {
			std::cerr << "Error: -k must be at most 6.\n";
			return 1;
		}




        std::cout << "Running OM-heap with w=" << w << " k=" << k << std::endl;
		std::cout << "Time limit (minutes): " << time_limit_minutes << std::endl;
		optimal_orders_routine(w, k, time_limit_minutes);
    }
    //else if (routine == "uhs_tests") {
    //    std::cout << "Running paper UHS tests" << std::endl;
    //    paper_optimal_uhs_routine();
    //}
    else if (routine == "calc_average_dens") {

		// Validate required args
		if (w == 0 || k == 0) {
			std::cerr << "Error: -w and -k are required parameters.\n"
				<< "Usage: ./my_program -routine <OM-heap|calc_average_dens> -w <w> -k <k> -time_limit_minutes <time_limit_in_minutes>\n";
			return 1;
		}
		// Validate required args
		if (w < 2) {
			std::cerr << "Error: -w must be at least 2.\n";
			return 1;
		}
		if (k < 1) {
			std::cerr << "Error: -k must be at least 1.\n";
			return 1;
		}
		if (k > 6) {
			std::cerr << "Error: -k must be at most 6.\n";
			return 1;
		}

        std::cout << "Running density calculation with w=" << w << " k=" << k << std::endl;
		random_density_routine(w, k);
    }
    else {
        std::cerr << "Unknown routine: " << routine
                  << "\nOptions: OM-heap | calc_average_dens\n";
        return 1;
    }

    return 0;
}






