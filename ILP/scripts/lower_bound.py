import numpy as np
import pandas as pd
import numpy as np
from decimal import Decimal, getcontext
getcontext().prec = 100  # You can increase this if needed


def find_first_mod_1_after_x(x, w):
    remainder = x % w
    if remainder == 0:
        to_add = 1
    else:
        to_add = w - remainder + 1
    return x + to_add


def is_prime(n):
    """Checks if a number is prime."""
    if n < 2:
        return False
    for i in range(2, int(Decimal(n).sqrt()) + 1):
        if n % i == 0:
            return False
    return True

def mobius(n):
    """Calculates the Möbius function value for n with high precision."""
    if n == 1:
        return 1

    # Use Decimal for high-precision sqrt calculation
    n_decimal = Decimal(n)
    sqrt_n = n_decimal.sqrt()

    p = 0  # Count of distinct prime factors


    # Iterate up to ceil(sqrt(n))
    for i in range(2,n):
        if n % i == 0:
            # If i is not prime, skip it
            if not is_prime(i):
                continue

            # If i divides n, check for square factors
            if (n // i) % i == 0:
                return 0  # n has a squared prime factor
            p += 1  # Count the distinct prime factor

    # Check if n itself is prime and greater than sqrt(n)
    if n > 1 and is_prime(n):
        p += 1

    # Return -1 if odd number of prime factors, 1 otherwise
    return -1 if p % 2 else 1

def aperiodic_necklaces(p, sigma):
    """
    Calculates the number of aperiodic necklaces of length p 
    over an alphabet of size sigma using the Möbius function.
    """
    return sum(mobius(p // d) * sigma**d for d in range(1, p + 1) if p % d == 0) // p

def get_lower_bound_improved(sigma, w, k):
    """Calculates the improved lower bound using high precision arithmetic."""
    # Convert inputs to native Python types
    w, k, sigma = int(w), int(k), int(sigma)

    # Use Decimal for high precision calculation of the first term
    first_term = Decimal(1) / (Decimal(sigma) ** Decimal(int(w + k)))

    # Get all divisors of (w + k)
    divisors = [i for i in range(1, w + k + 1) if (w + k) % i == 0]

    
    sum_terms = sum(aperiodic_necklaces(d, sigma) * np.ceil(d / w) for d in divisors)

    return float(first_term * Decimal(sum_terms))


def get_lower_bound_improved_number_windows(sigma, w, k):
    """Calculates the improved lower bound using high precision arithmetic."""
    # Convert inputs to native Python types
    w, k, sigma = int(w), int(k), int(sigma)


    # Get all divisors of (w + k)
    divisors = [i for i in range(1, w + k + 1) if (w + k) % i == 0]

    
    sum_terms = sum(aperiodic_necklaces(d, sigma) * np.ceil(d / w) for d in divisors)

    return sum_terms



def get_best_lower_bound_number_windows(sigma,w,k):
    bound_1 = get_lower_bound_improved_number_windows(sigma,w,k)
    #print(f"Lower bound option 1: {bound_1}")
    k_2 = find_first_mod_1_after_x(k,w)
    k_diff = k_2 - k
    bound_2 = get_lower_bound_improved_number_windows(sigma,w,k_2) / (sigma ** k_diff)
    #print(f"Lower bound option 2 (k adjusted to {k_2}): {bound_2}")
    return max(bound_1, bound_2)