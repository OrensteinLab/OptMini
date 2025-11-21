from pyexpat import model
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import json
import itertools
from scripts.scripts import *
from scripts.load_gm_order import *


def get_model():
    with open('gurobi.json', 'r') as json_file:
        params = json.load(json_file)
        
        env = gp.Env(params=params)

        # Create the model within the Gurobi environment
        model = gp.Model('Minimizer ILP', env=env)

        return model
    

def set_constraints(model, sigma, w, k, use_gm_as_starting_point = True, extra_constraints = None):
    epsilon = 1  # was 1e-3




    
    # Generate all possible k-mers
    kmer_list = [''.join(p) for p in itertools.product(map(str, range(sigma)), repeat=k)]
    num_kmers = len(kmer_list)

    # Map k-mers to indices for easy lookup
    kmer_to_index = {kmer: idx for idx, kmer in enumerate(kmer_list)}

    # Generate all possible sequences of length k + w
    sequence_length = k + w
    sequence_list = [''.join(p) for p in itertools.product(map(str, range(sigma)), repeat=sequence_length)]
    num_sequences = len(sequence_list)

    n_kmers=  sigma**k
    #spread_factor = 1000 # spread factor 1 is default
    #new_ub = spread_factor * (n_kmers)
    new_ub = n_kmers

    x_vars = model.addVars(num_kmers, vtype=GRB.INTEGER, lb=0, ub=new_ub, name="x")
    y_vars = model.addVars(num_sequences, vtype=GRB.BINARY, name="y")



    if does_gm_order_exist(w,k,sigma) and use_gm_as_starting_point:
        print("Using GreedyMini order as starting point")
        loaded_x = load_gm_order(w,k, sigma)
        # Set the initial values for the variables using the `start` attribute
        model.update()  # Ensure all variables are added to the model before setting the start values
        for i in range(len(loaded_x)):
            #x_vars[i].start = loaded_x[i]* spread_factor  # Initialize the search with the loaded values
            x_vars[i].start = loaded_x[i]
    else:
        print("Using random starting point")




    


    # Set the objective: Minimize sum of y_u
    model.setObjective(
        gp.quicksum(y_vars[u] for u in range(num_sequences)),
        GRB.MINIMIZE
    )

    # Iterate over all sequences to define constraints
    for u_idx, u in enumerate(sequence_list):
        # Extract k-mers in the sequence
        kmers_in_u = [u[i:i+k] for i in range(len(u) - k + 1)]
        unique_kmers = list(dict.fromkeys(kmers_in_u))  # Remove duplicates, keep order

        # Indices of all k-mers in the sequence
        kmer_indices = [kmer_to_index[kmer] for kmer in kmers_in_u]
        unique_kmer_indices = [kmer_to_index[kmer] for kmer in unique_kmers]

        # Indices of the first and last k-mers
        first_kmer_idx = kmer_to_index[kmers_in_u[0]]
        last_kmer_idx = kmer_to_index[kmers_in_u[-1]]

        # Check if the last k-mer appears more than once
        last_appears_twice = kmer_indices.count(last_kmer_idx) > 1

        # **Singleton Sequence Handling**
        if len(unique_kmers) == 1:
            # If there's only one unique k-mer, y_u must be 1
            model.addConstr(y_vars[u_idx] == 1, name=f"singleton_yu{u_idx}")
            continue  # Skip to the next sequence



        ### **Condition 1: First k-mer is greater than at least one other unique k-mer**
        other_kmer_indices_first = [k_idx for k_idx in unique_kmer_indices if k_idx != first_kmer_idx]

        if other_kmer_indices_first:
            # Compute the minimum rank among the other k-mers
            min_other_first = model.addVar(lb=0, ub=new_ub, name=f"min_first_others_{u_idx}")
            model.addGenConstrMin(min_other_first, [x_vars[k] for k in other_kmer_indices_first])

            # Create binary z_first_u = 1 if min_other <= first_kmer - epsilon
            z_first_u = model.addVar(vtype=GRB.BINARY, name=f"z_first_u{u_idx}")
            model.addGenConstrIndicator(
                z_first_u, True,
                min_other_first <= x_vars[first_kmer_idx] - epsilon,
                name=f"z1_first_le_min_{u_idx}"
            )
            model.addGenConstrIndicator(
                z_first_u, False,
                min_other_first >= x_vars[first_kmer_idx],
                name=f"z0_first_ge_min_{u_idx}"
            )
        else:
            # Only one unique k-mer — skip
            z_first_u = model.addVar(vtype=GRB.BINARY, name=f"z_first_u{u_idx}")
            model.addConstr(z_first_u == 0, name=f"z_first_singleton_{u_idx}")



        ### **Linking Conditions to y_u**

        if last_appears_twice:
            # If the suffix isn't unique, only consider the first condition
            s_u = z_first_u
        else:
            ### **Condition 2: Last k-mer is greater than at least one other unique k-mer**
            other_kmer_indices_last = [k_idx for k_idx in unique_kmer_indices if k_idx != last_kmer_idx]

            if other_kmer_indices_last:
                # Compute the minimum rank among all non-last k-mers
                min_other_last = model.addVar(lb=0, ub=new_ub, name=f"min_last_others_{u_idx}")
                model.addGenConstrMin(min_other_last, [x_vars[k] for k in other_kmer_indices_last])

                # Binary z_last_u = 1 if min_other_last <= x_vars[last_kmer_idx] - epsilon
                z_last_u = model.addVar(vtype=GRB.BINARY, name=f"z_last_u{u_idx}")
                model.addGenConstrIndicator(
                    z_last_u, True,
                    min_other_last <= x_vars[last_kmer_idx] - epsilon,
                    name=f"z1_last_le_min_{u_idx}"
                )
                model.addGenConstrIndicator(
                    z_last_u, False,
                    min_other_last >= x_vars[last_kmer_idx],
                    name=f"z0_last_ge_min_{u_idx}"
                )
            else:
                # Only one unique k-mer (or no valid others)
                z_last_u = model.addVar(vtype=GRB.BINARY, name=f"z_last_u{u_idx}")
                model.addConstr(z_last_u == 0, name=f"z_last_singleton_{u_idx}")



            # Combine conditions using logical AND
            s_u = model.addVar(vtype=GRB.BINARY, name=f"s_u{u_idx}")
            model.addGenConstrAnd(s_u, [z_first_u, z_last_u], name=f"s_u_and_{u_idx}")


        # Set y_u >= 1 - s_u
        model.addConstr(y_vars[u_idx] >= 1 - s_u,
                        name=f"y_constraint_u{u_idx}")
        


        


def model_optimize(model, lower_bound, thread_count = 8, time_limit_in_mins = 0, seed = 42):

    model.setParam('Threads', thread_count) 
    model.setParam('Seed', seed) 

    if time_limit_in_mins > 0:
        model.setParam(GRB.Param.TimeLimit, time_limit_in_mins*60)


    # Define the callback function to stop early if the lower bound is reached
    def callback(model, where):
        if where == GRB.Callback.MIP:
            # Get the best objective value found so far
            obj_value = model.cbGet(GRB.Callback.MIP_OBJBST)
            if obj_value <= lower_bound:
                print(f"Early stop: Objective value {obj_value} reached lower bound {lower_bound}.")
                model.terminate()


    model.optimize(callback)


    
    return model
    

def validate_order(order):
    # Step 1: Pair each element with its index
    value_index_list = list(enumerate(order))

    # Step 2: Sort the list based on values (ties broken by index)
    sorted_value_index_list = sorted(value_index_list, key=lambda x: x[1])

    # Step 3: Assign ranks and place them back into a list matching the original order
    ranks = [0] * len(order)
    for rank, (index, value) in enumerate(sorted_value_index_list, start=1):
        ranks[index] = rank

    return ranks



def extract_order(model):
    x_values = []
    for var in model.getVars():
        if var.varName[0] == 'x':
            x_values.append(var.x)

    #print("Optimal x values:")
    #print(x_values)
    return x_values


def get_density(model,sigma, w, k):
    sum_y = 0
    for var in model.getVars():
        if var.varName[0] == 'y':
            if var.x > 0.5:
                sum_y += 1
    
    density = sum_y / (sigma**(w+k))
    print(f"Density: {density}")

    density_factor = (1+w) * density
    print(f"Density factor: {density_factor}")

    return density



def get_density_new_method(model, sigma, w,k):
    return get_density(model, sigma, w, k)



def get_density_old_method(model, sigma, w,k):
    order = extract_order(model)
    new_order = validate_order(order)
    charge_count = prob_gc(w,k,new_order)
    density = charge_count / (sigma**(w+k))
    print(f"Density: {density}")

    density_factor = (1+w) * density
    print(f"Density factor: {density_factor}")

    return density
