
#include "functions.h" 







std::vector<uint64_t> convert_to_order(std::vector<uint64_t> original_order, uint32_t k) {
	uint64_t n_kmers = static_cast<uint64_t>(1) << k;
	uint64_t default_value = static_cast<uint64_t>(1) << k;
	std::vector<uint64_t> new_order = std::vector<uint64_t>(n_kmers, default_value);

	for (uint64_t i = 0; i < original_order.size(); i++) {
		new_order[original_order[i]] = i;
	}
	return new_order;
}



//std::vector<bool> get_bits(uint64_t n_bits) {
//	std::vector<bool> bits(n_bits);
//
//	std::random_device rd;
//	std::mt19937 gen(rd());
//	std::uniform_int_distribution<int> dist(0, 1);
//
//	for (size_t i = 0; i < n_bits; ++i) {
//		bits[i] = dist(gen); // Either 0 or 1
//	}
//
//	return bits;
//
//}
