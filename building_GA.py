import numpy as np
import math,random
# from encode import encode_float_vector_to_bitstream
# from decode import decode_bitstream_to_float_vector
# def encode_float_vector_to_bitstream(float_vector, search_space):
#     normalized_vector = (float_vector - search_space[:, 0]) / (search_space[:, 1] - search_space[:, 0])
#     bitstream = np.round(normalized_vector).astype(int)
#     return bitstream

# def decode_bitstream_to_float_vector(bitstream, search_space):
#     return search_space[:, 0] + bitstream * (search_space[:, 1] - search_space[:, 0])
# def encode_float_vector_to_bitstream(float_vector, search_space, num_genes):
#     bitstream = np.zeros(num_genes, dtype=int)
#     for i in range(len(search_space)):
#         lower_bound, upper_bound = search_space[i]
#         gene_range = upper_bound - lower_bound
#         gene_value = (float_vector[i] - lower_bound) / gene_range
#         gene_index = int(np.round(gene_value * (num_genes - 1)))
#         bitstream[i] = gene_index
#     return bitstream

# def decode_bitstream_to_float_vector(bitstream, search_space, num_genes):
#     float_vector = np.zeros(len(search_space), dtype=float)
#     for i in range(len(search_space)):
#         lower_bound, upper_bound = search_space[i]
#         gene_range = upper_bound - lower_bound
#         gene_value = bitstream[i] / (num_genes - 1)
#         float_vector[i] = lower_bound + gene_value * gene_range
#     return float_vector

import numpy as np

def float_vector_to_binary(float_vector, search_space):
    binary_vector = ""
    
    for i in range(len(float_vector)):
        # Get the range for the current gene
        lower_bound, upper_bound = search_space[i]
        
        # Calculate the range of the current gene
        gene_range = upper_bound - lower_bound
        
        # Calculate the number of bits needed to represent the gene value
        num_bits = int(np.ceil(np.log2(gene_range + 1)))
        
        # Normalize the gene value within [0, 1]
        normalized_value = (float_vector[i] - lower_bound) / gene_range
        
        # Convert the normalized value to binary with the required number of bits
        binary_representation = format(int(normalized_value * (2 ** num_bits - 1)), f"0{num_bits}b")
        
        # Append the binary representation of the current gene to the binary vector
        binary_vector += binary_representation
    
    return binary_vector

# Example usage
float_vector = [10.5, 15.0, 6.8]  # Replace with your float vector
search_space = [(8.0, 20.0), (8.0, 20.0), (6.0, 18.0)]  # Replace with your search space
binary_representation = float_vector_to_binary(float_vector, search_space)
print(binary_representation)


def encode_float_vector_to_bitstream(float_vector, search_space, num_genes): #encode gets the initial float vector containing the genes, number of genes (equal to length variables and bounds for each variable
    bitstream = np.zeros(num_genes, dtype=int) #bitstream is the number of genes
    for i in range(len(search_space)):                 
        lower_bound, upper_bound = search_space[i]  #define upper and lower bound of genes
        if lower_bound == upper_bound:         
            gene_value = 0.5  # Assign a default value when the range is zero
        else:
            gene_range = upper_bound - lower_bound
            gene_value = (float_vector[i] - lower_bound) / gene_range  #normalize genes according gene range
        gene_index = int(np.round(gene_value * (num_genes - 1)))
        bitstream[i] = gene_index
    return bitstream


def decode_bitstream_to_float_vector(bitstream, search_space, num_genes):
    float_vector = np.zeros(len(search_space), dtype=float)
    for i in range(len(search_space)):
        lower_bound, upper_bound = search_space[i]
        gene_range = upper_bound - lower_bound
        gene_value = bitstream[i] / (num_genes - 1)
        float_vector[i] = lower_bound + gene_value * gene_range
    return float_vector

def initialize_population(population_size, search_space, num_genes):
    population = np.empty((population_size, num_genes), dtype=int)
    for i in range(population_size):
        random_float_vector = np.random.uniform(search_space[:, 0], search_space[:, 1])
        encoded_individual = encode_float_vector_to_bitstream(random_float_vector, search_space, num_genes)
        population[i] = encoded_individual
    return population


def calculate_num_bits_per_element(search_space, num_genes):
    num_bits_per_element = []
    for i in range(num_genes):
        lower_bound, upper_bound = search_space[i]
        gene_range = upper_bound - lower_bound
        num_bits = int(math.ceil(math.log2(num_genes)))
        num_bits_per_element.append(num_bits)
    return num_bits_per_element

# def initialize_population(population_size, search_space):
#     return np.random.randint(2, size=(population_size, len(search_space)))

# def initialize_population(population_size, num_genes):
    # return np.random.randint(2, size=(population_size, num_genes))

def apply_constraints(float_solution, search_space):
    constrained_solution = np.clip(float_solution, search_space[:, 0], search_space[:, 1])
    return constrained_solution

def objective_function(encoded_solution, search_space):
    float_solution = decode_bitstream_to_float_vector(encoded_solution, search_space, num_genes)
    constrained_solution = apply_constraints(float_solution, search_space)

    XL, XB1, XB2, XB3, XBT, XK, XH = constrained_solution[0], constrained_solution[2], constrained_solution[4], constrained_solution[6], constrained_solution[8], constrained_solution[10], constrained_solution[12]
    YL, YB1, YB2, YB3, YBT, YK, YH = constrained_solution[1], constrained_solution[3], constrained_solution[5], constrained_solution[7], constrained_solution[9], constrained_solution[11], constrained_solution[13]
    penalty=0
    # Your objective function here
    F = XL * YL + XB1 * YB1 + XB2 * YB2 + XB3 * YB3 + 2 * XBT * YBT + 2 * XK * YK + XH * YH + 3
    if not (120 <= XL * YL <= 300):
        penalty += abs(min(120 - XL * YL, XL * YL - 300))
    if not (XL * YL / F == 1.5):
        penalty += abs(XL * YL / F - 1.5)
# ------------------------------------------------------------------------------------------------------

    if not (50 <= XK * YK <= 120):
        penalty += abs(min(50 - XK * YK, XK * YK - 120))
# -----------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------

    if not (19 <= XH * YH <= 72):
        penalty += abs(min(19 - XH * YH, XH * YH - 72))
    # ------------------------------------------------------------------------------------------------------

    if not (100 <= XB1 * YB1 <= 180):
        penalty += abs(min(100 - XB1 * YB1, XB1 * YB1 - 180))
    if not (XB1 * YB1 / F == 1.5):
        penalty += abs(XB1 * YB1 / F - 1.5)
    # ---------------------------------------------------------------------------------

    if not (100 <= XB2 * YB2 <= 180):
        penalty += abs(min(100 - XB2 * YB2, XB2 * YB2 - 180))
    if not (XB1 * YB1 / F == 1.5):
        penalty += abs(XB2 * YB2 / F - 1.5)
    # -----------------------------------------------------------------------------------

    if not (100 <= XB3 * YB3 <= 180):
        penalty += abs(min(100 - XB3 * YB3, XB3 * YB3 - 180))
    if not (XB3 * YB3 / F == 1.5):
        penalty += abs(XB3 * YB3 / F - 1.5)
    F+=penalty

    return F

def crossover(parent1, parent2, crossover_rate):
    if np.random.rand() < crossover_rate:
        crossover_point = np.random.randint(1, len(parent1) - 1)
        child = np.concatenate((parent1[:crossover_point], parent2[crossover_point:]))
    else:
        child = parent1.copy()
    return child

def mutate(solution, mutation_rate):
    mutated_solution = solution.copy()
    mutation_mask = np.random.rand(len(solution)) < mutation_rate
    mutated_solution[mutation_mask] = 1 - mutated_solution[mutation_mask]
    return mutated_solution
def genetic_algorithm(population_size, generations, mutation_rate, crossover_rate, search_space, num_genes):
    population = initialize_population(population_size, search_space, num_genes)

    for generation in range(generations):
        # Evaluate the fitness of each individual
        fitness_values = np.array([-objective_function(decode_bitstream_to_float_vector(individual, search_space, num_genes), search_space) for individual in population])

        # Select parents based on fitness
        parents_indices = np.argsort(fitness_values)[:2]
        parent1, parent2 = population[parents_indices]

        # Crossover to create offspring
        child = crossover(parent1, parent2, crossover_rate)

        # Mutation
        child = mutate(child, mutation_rate)

        # Replace the least fit individual with the offspring
        population[np.argmax(fitness_values)] = child

        # Print the best solution in each generation
        best_solution = population[np.argmin(fitness_values)]
        best_fitness = objective_function(decode_bitstream_to_float_vector(best_solution, search_space, num_genes), search_space)
        print(f"Generation {generation + 1}, Best Fitness: {best_fitness}")

    # Return the best solution found
    best_solution = population[np.argmin(fitness_values)]
    return best_solution

# Define the search space
search_space = np.array([[8, 20], [8, 20], [10, 17], [10, 17], [9, 20], [9, 20], [8, 18], [8, 18], [5.5, 5.5], [8.5, 8.5], [6, 18], [6, 18], [5.5, 5.5], [3.5, 6]])
# Set algorithm parameters
a=3
if a==1:
    population_size = 50 #large population size
    generations = 50
    mutation_rate = 0.009
    crossover_rate = 0.85
    num_genes = len(search_space) # Number of genes per individual
if a==2:
    population_size = 100 #large population size
    generations = 50
    mutation_rate = 0.006
    crossover_rate = 0.90
    num_genes = len(search_space) # Number of genes per individual
if a==3:
    population_size = 200 #large population size
    generations = 50
    mutation_rate = 0.002
    crossover_rate = 0.95
    num_genes = len(search_space) # Number of genes per individual


# Run the genetic algorithm
best_solution = genetic_algorithm(population_size, generations, mutation_rate, crossover_rate, search_space, num_genes)
# Calculate the number of bits required for each element
num_bits_per_element = calculate_num_bits_per_element(search_space, num_genes)
# Calculate the total number of bits for the entire bitstream
num_bits_total = sum(num_bits_per_element)

# Decode the best solution to get the final result
final_result = apply_constraints(decode_bitstream_to_float_vector(best_solution, search_space, num_genes), search_space)
print("Final Result:", final_result)
print('number of bits in each element should be{}, and the percision is{}'.format(num_bits_total,1/num_bits_total))
obj_val=objective_function(best_solution,search_space)
print('the objective value is:{}'.format(obj_val))
a=[random.uniform(_[0],_[1]) for _ in search_space]
b=encode_float_vector_to_bitstream(a,search_space,num_genes)
binary_bitstream = [bin(value)[2:].zfill(4) for value in b]

# c=decode_bitstream_to_float_vector()
print('the first encoded value is{} and the decoded value is{}'.format(binary_bitstream,a) )
a=[random.uniform(_[0],_[1]) for _ in search_space]
b=encode_float_vector_to_bitstream(a,search_space,num_genes)
binary_bitstream = [bin(value)[2:].zfill(4) for value in b]
print('the first encoded value is{} and the decoded value is{}'.format(binary_bitstream,a) )
a=best_solution
binary_bitstream = [bin(value)[2:].zfill(4) for value in b]
print('the first encoded value is{} and the decoded value is{}'.format(binary_bitstream,decode_bitstream_to_float_vector(a,search_space,num_genes)) )




#----------performance-----++++++++++++++++++------------------------

import numpy as np

# Define your GA functions here (encode, decode, objective_function, etc.)

# Set the number of runs to assess reliability
num_runs = 1
results = []

for run in range(num_runs):
    # Run the GA
    best_solution = genetic_algorithm(population_size, generations, mutation_rate, crossover_rate, search_space, num_genes)
    
    # Calculate objective value and other statistics
    obj_val = objective_function(decode_bitstream_to_float_vector(best_solution, search_space, num_genes), search_space)
    
    # Store the results
    results.append({
        "Run": run + 1,
        "Objective Value": obj_val,
        # Add more statistics or measurements here
    })

# Calculate statistics across runs
objective_values = [result["Objective Value"] for result in results]
avg_obj_value = np.mean(objective_values)
std_dev_obj_value = np.std(objective_values)

# Measure time for each run (optional)
import time
runtimes = []

for run in range(num_runs):
    start_time = time.time()
    best_solution = genetic_algorithm(population_size, generations, mutation_rate, crossover_rate, search_space, num_genes)
    end_time = time.time()
    runtime = end_time - start_time
    runtimes.append(runtime)

avg_runtime = np.mean(runtimes)
std_dev_runtime = np.std(runtimes)

# Print and analyze the results
print("Results for", num_runs, "runs of the GA:")
print("Average Objective Value:", avg_obj_value)
print("Standard Deviation of Objective Value:", std_dev_obj_value)
print("Average Runtime (seconds):", avg_runtime)
print("Standard Deviation of Runtime (seconds):", std_dev_runtime)

# You can add more analysis here to assess reliability and other factors.



# ---------------plot the performance-----++++------------------------

import numpy as np
import matplotlib.pyplot as plt

# Define your GA functions here (encode, decode, objective_function, etc.)

# Set the number of runs to assess reliability
num_runs = 1
results = []

for run in range(num_runs):
    # Run the GA
    best_solution = genetic_algorithm(population_size, generations, mutation_rate, crossover_rate, search_space, num_genes)
    
    # Calculate objective value and other statistics
    obj_val = objective_function(decode_bitstream_to_float_vector(best_solution, search_space, num_genes), search_space)
    
    # Store the results
    results.append({
        "Run": run + 1,
        "Objective Value": obj_val,
        # Add more statistics or measurements here
    })

# Extract objective values and runtimes
objective_values = [result["Objective Value"] for result in results]
runtimes = []

for run in range(num_runs):
    start_time = time.time()
    best_solution = genetic_algorithm(population_size, generations, mutation_rate, crossover_rate, search_space, num_genes)
    end_time = time.time()
    runtime = end_time - start_time
    runtimes.append(runtime)

# Create a performance plot for objective values
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(range(1, num_runs + 1), objective_values, marker='o', linestyle='-')
plt.title("Objective Value vs. Run")
plt.xlabel("Run")
plt.ylabel("Objective Value")

# Create a performance plot for runtimes
plt.subplot(1, 2, 2)
plt.plot(range(1, num_runs + 1), runtimes, marker='o', linestyle='-')
plt.title("Runtime vs. Run")
plt.xlabel("Run")
plt.ylabel("Runtime (seconds)")

plt.tight_layout()
plt.show()
