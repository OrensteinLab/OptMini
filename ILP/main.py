import argparse
from scripts.ilp import *
from scripts.lower_bound import *


def main(w, k, sigma, time_limit, threads):
    lower_bound = get_best_lower_bound_number_windows(sigma, w, k)
    print(f"Lower bound: {lower_bound}")


    model = get_model()
    set_constraints(model, sigma, w, k)

    # Check if tuning found a set of parameters
    # if model.tuneResultCount > 0:
    #     # Load the best set of parameters
    #     model.getTuneResult(0)

    model.setParam('Presolve', 2)  # Enable aggressive presolve
    #model.setParam('Cuts', 2)  # Enable aggressive cut generation

    model_optimize(model,lower_bound, threads, time_limit)
    density_2 = get_density_new_method(model, sigma, w, k)

    #order = extract_order(model, sigma, w, k)
    order = extract_order(model)

    # save order to file (it's a list)
    # TODO: only do this if optimal solution is found
    with open(f"order_w{w}_k{k}_sigma{sigma}.txt", "w") as f:
        f.write(str(order))
        print(f"Order saved to order_w{w}_k{k}_sigma{sigma}.txt")



    if sigma == 2:
        density_1 = get_density_old_method(model, sigma, w, k)
        print(f"Density calculated with old method: {density_1}")
        print(f"Density calculated with ILP variables: {density_2}")
        return density_1, density_2
    else:
        return None, density_2
    
    

import argparse
import os
import multiprocessing

if __name__ == '__main__':
    # Determine the number of physical cores
    default_threads = multiprocessing.cpu_count() // 2

    # Create the parser
    parser = argparse.ArgumentParser(description="Process input arguments for the program.")

    # Add arguments
    parser.add_argument('-w', type=int, required=True, help='Window size (w)')
    parser.add_argument('-k', type=int, required=True, help='K-mer length (k)')
    parser.add_argument('-sigma', type=int, default=2, help='Sigma value (default: 2)')
    parser.add_argument('-time_limit', type=int, default=1, help='Time limit in minutes (default: 1)')
    parser.add_argument('-threads', type=int, default=default_threads, help=f'Number of threads to use (default: {default_threads})')

    # Parse the arguments
    args = parser.parse_args()

    # Call main with parsed arguments
    main(args.w, args.k, args.sigma, args.time_limit, args.threads)
