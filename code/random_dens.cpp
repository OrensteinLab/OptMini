
#include "random_dens.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
using boost::multiprecision::cpp_dec_float_100;
using boost::multiprecision::uint256_t;

std::vector<double> randens_w(uint32_t k, uint32_t w) {
	const uint64_t n_kmers = 1ULL << k;
	const uint64_t n_masks = 1ULL << n_kmers;



	// Transition table
	std::vector<std::array<uint32_t, 2>> table(n_kmers);

	std::vector<std::vector<uint64_t>> tot_gc(w, std::vector<uint64_t>(n_kmers, 0ULL));

	std::vector<uint64_t> denom(n_kmers, 0ULL);
	denom[0] = n_kmers;
	for (uint64_t i = 1; i < n_kmers; i++)
		denom[i] = denom[i - 1] * (n_kmers - i) / i;

	for (uint64_t mask = 0; mask < n_masks; ++mask) {
		uint32_t rank = std::popcount(mask); // bit count
		uint32_t c0 = rank, c1 = 0;
		std::vector<uint32_t> order(n_kmers);

		for (uint64_t kmer = 0; kmer < n_kmers; ++kmer)
			if ((mask >> kmer) & 1ULL) {
				order[kmer] = c1;
				c1++;
			}
			else {
				order[kmer] = c0;
				c0++;
			}
			

		for (uint64_t i = 0; i < n_kmers; ++i) {
			uint64_t a = (i << 1) & (n_kmers - 1);
			table[order[i]][0] = order[a];
			table[order[i]][1] = order[a + 1];
		}

		for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
			if (((mask >> kmer) & 1ULL) == 0) {
				std::vector<uint64_t> old(n_kmers, 0), New(n_kmers, 0);
				uint64_t kmer_rank = order[kmer];
				old[table[kmer_rank][0]] = 1;
				old[table[kmer_rank][1]] = 1;

				for (uint32_t i = 1; i < w; ++i) { // DP for prefix gc
					for (uint64_t ii = rank; ii < n_kmers; ++ii) {
						New[table[ii][0]] += old[ii];
						New[table[ii][1]] += old[ii];
						old[ii] = 0;
					}
					std::swap(old, New);
					uint64_t g = 0;
					for (uint64_t ii = rank; ii < n_kmers; ++ii) g += old[ii];
					if (g == 0) break;
					tot_gc[i][rank] += g;
				}

				std::fill(old.begin(), old.end(), 1);
				std::fill(New.begin(), New.end(), 0);

				for (uint32_t i = 0; i < w; ++i) {
					old[order[kmer]] = 0;
					for (uint64_t ii = rank; ii < n_kmers; ++ii) {
						New[table[ii][0]] += old[ii];
						New[table[ii][1]] += old[ii];
						old[ii] = 0;
					}
					std::swap(old, New);
					uint64_t g = 0;
					for (uint64_t ii = rank; ii < n_kmers; ++ii) g += old[ii];
					if (g == 0) break;
					tot_gc[i][rank] += old[order[kmer]];
				}
			}
		}
	}

	// ====== High-precision density calculation ======
	std::vector<double> ans(w - 1, 0.0);
	for (uint32_t i = 1; i < w; ++i) {
		cpp_dec_float_100 temp = 0;
		for (uint64_t rank = 0; rank < n_kmers; ++rank)
			temp += cpp_dec_float_100(tot_gc[i][rank]) / cpp_dec_float_100(denom[rank]);

		cpp_dec_float_100 density = temp / pow(cpp_dec_float_100(2), k + i + 1) * cpp_dec_float_100(i + 2);
		ans[i - 1] = static_cast<double>(density);
	}

	return ans;
}



// same but uses 256 bits for operations 
std::vector<double> randens_w_256(uint32_t k, uint32_t w) {
	const uint64_t n_kmers = 1ULL << k;
	const uint64_t n_masks = 1ULL << n_kmers;


	// compute final mask value (inclusive)
	uint64_t n_masks_minus_one = 0;
	if (k < 6)
		n_masks_minus_one = (1ULL << n_kmers) - 1ULL;
	else if (k == 6)
		n_masks_minus_one = std::numeric_limits<uint64_t>::max();
	else
		throw std::invalid_argument("k must be <= 6");

	// Transition table
	std::vector<std::array<uint32_t, 2>> table(n_kmers);

	std::vector<std::vector<uint256_t>> tot_gc(w, std::vector<uint256_t>(n_kmers, 0ULL));

	std::vector<uint256_t> denom(n_kmers, 0ULL);
	denom[0] = n_kmers;
	for (uint64_t i = 1; i < n_kmers; i++)
		denom[i] = denom[i - 1] * (n_kmers - i) / i;

	for (uint64_t mask = 0; mask<n_masks; ++mask) {
		uint32_t rank = std::popcount(mask);
		uint32_t c0 = rank, c1 = 0;
		std::vector<uint32_t> order(n_kmers);

		for (uint64_t kmer = 0; kmer < n_kmers; ++kmer)
			order[kmer] = ((mask >> kmer) & 1ULL) ? c1++ : c0++;

		for (uint64_t i = 0; i < n_kmers; ++i) {
			uint64_t a = (i << 1) & (n_kmers - 1);
			table[order[i]][0] = order[a];
			table[order[i]][1] = order[a + 1];
		}

		for (uint64_t kmer = 0; kmer < n_kmers; ++kmer) {
			if (((mask >> kmer) & 1ULL) == 0) {
				std::vector<uint256_t> old(n_kmers, 0), New(n_kmers, 0);

				old[table[order[kmer]][0]] = 1;
				old[table[order[kmer]][1]] = 1;

				for (uint32_t i = 1; i < w; ++i) {
					for (uint64_t ii = rank; ii < n_kmers; ++ii) {
						New[table[ii][0]] += old[ii];
						New[table[ii][1]] += old[ii];
						old[ii] = 0;
					}
					std::swap(old, New);
					uint256_t g = 0;
					for (uint64_t ii = rank; ii < n_kmers; ++ii)
					{
						g += old[ii];
					}
					if (g == 0)
					{
						break;
					}
					tot_gc[i][rank] += g;
				}

				std::fill(old.begin(), old.end(), 1);
				std::fill(New.begin(), New.end(), 0);

				for (uint32_t i = 0; i < w; ++i) {
					old[order[kmer]] = 0;
					for (uint64_t ii = rank; ii < n_kmers; ++ii) {
						New[table[ii][0]] += old[ii];
						New[table[ii][1]] += old[ii];
						old[ii] = 0;
					}
					std::swap(old, New);
					uint256_t g = 0;
					for (uint64_t ii = rank; ii < n_kmers; ++ii) {
						g += old[ii];
					}
					if (g == 0)
					{
						break;
					}
					tot_gc[i][rank] += old[order[kmer]];
				}
			}
		}
	}

	// ====== High-precision density calculation ======
	std::vector<double> ans(w - 1, 0.0);
	for (uint32_t i = 1; i < w; ++i) {
		cpp_dec_float_100 temp = 0;
		for (uint64_t rank = 0; rank < n_kmers; ++rank)
			temp += cpp_dec_float_100(tot_gc[i][rank]) / cpp_dec_float_100(denom[rank]);

		cpp_dec_float_100 density = temp / pow(cpp_dec_float_100(2), k + i + 1) * cpp_dec_float_100(i + 2);
		ans[i - 1] = static_cast<double>(density);
	}

	return ans;
}




double get_random_density(uint32_t w, uint32_t k) {
	std::vector<double> densities = randens_w(k, w);
	// print the size of densities
	std::cout << "Size of densities: " << densities.size() << std::endl;
	// print w-2
	std::cout << "w-2: " << (w - 2) << std::endl;
	return randens_w(k, w)[w - 2];
}



double get_random_density_256(uint32_t w, uint32_t k) {
	return randens_w_256(k, w)[w - 2];
}
