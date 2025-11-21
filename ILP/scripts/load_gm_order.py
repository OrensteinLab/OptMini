import struct







def get_explicitly_extended_order(binary_order):
    original_size = len(binary_order)
    new_size =  original_size ** 2 #from binary ABV to DNA
    dna_order = [0] * new_size


    # Loop through each index in the new order
    for i in range(new_size):
        # Convert the index to a binary string, padded with leading zeros
        bit_string = f"{i:0{new_size.bit_length()-1}b}"

        # Extract odd and even bits from the bit string
        odd_bits = bit_string[::2]   # Bits at odd positions (1-based index)
        even_bits = bit_string[1::2]   # Bits at even positions

        # Reconstruct numbers from the odd and even bits
        odd_number = int(odd_bits, 2) if odd_bits else 0
        even_number = int(even_bits, 2) if even_bits else 0

        dna_order[i] = binary_order[odd_number] * original_size + even_number

        if dna_order[i] > new_size:
            dna_order[i] = new_size
    return dna_order



def does_gm_order_exist(w, k, sigma):
    filename = f'greedymini orders/w{w}_k{k}.bin'
    if sigma in [2, 4]:
        try:
            with open(filename, 'rb'):
                return True
        except FileNotFoundError:
            return False


def load_gm_order(w,k, sigma):
    filename = f'greedymini orders/w{w}_k{k}.bin'

    if sigma == 2:
        order = load_vector_from_file(filename)
    elif sigma == 4:
        order = load_vector_from_file(filename)
        order = get_explicitly_extended_order(order)
    else:
        raise ValueError(f"Invalid sigma value: {sigma}")


    return order

def load_vector_from_file(filename):
    with open(filename, 'rb') as file:
        # Read the size of the vector (stored as size_t in C++, which is platform-dependent)
        size_data = file.read(struct.calcsize('P'))  
        size = struct.unpack('P', size_data)[0]
        
        # Read the rest of the data as uint64_t values
        vector_data = file.read(size * struct.calcsize('Q')) 
        vector = struct.unpack(f'{size}Q', vector_data)
        
        return list(vector)

