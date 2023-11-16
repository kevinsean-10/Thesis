from random import choices, randint, randrange, random
from typing import List, Optional, Callable, Tuple
import numpy as np
import pandas as pd
import sobol_seq
import matplotlib.pyplot as plt
from functools import partial

Genome = List[int]
Population = List[Genome]
FitnessFunc = Callable[[Genome], float]
PopulateFunc = Callable[[], Population]
SelectionFunc = Callable[[Population, FitnessFunc], Tuple[Genome, Genome]]
CrossoverFunc = Callable[[Genome, Genome], Tuple[Genome, Genome]]
MutationFunc = Callable[[Genome], Genome]
PrinterFunc = Callable[[Population, int, FitnessFunc], None]

"""GENERATE POINTS USING SOBOL SEQUENCE"""
def generate_points(dim,npoint,low=-10,high=10):
    if type(low) != type(high):
        raise TypeError('The type of "low" and "high" should be the same.')
    if type(low) == int:
        boundaries = [(low,high) for _ in range (dim)]
    elif type(low) == list or type(low) == np.ndarray:
        if len(low) != len(high):
            raise TypeError('The length of "low" and "high" should be the same.')
        else:
            boundaries = [(low[i],high[i]) for i in range (len(low))]

    # Generate Sobol sequence points
    sobol_points = sobol_seq.i4_sobol_generate(dim, npoint)

    # Scale the Sobol points to fit within the specified boundaries
    scaled_points = []
    for i in range(dim):
        a, b = boundaries[i]
        scaled_dim = a + sobol_points[:, i] * (b - a)
        scaled_points.append(scaled_dim)

    # Transpose the scaled points to get points per dimension
    scaled_points = np.array(list(map(list, zip(*scaled_points))))
    return scaled_points

def encode_number(number, min_value=-10, max_value=10, num_bits=32):
    # Normalize the number to a value between 0 and 1
    normalized = (number - min_value) / (max_value - min_value)
    # Convert it to an integer representation
    int_representation = int(normalized * (2**num_bits - 1))
    # Convert the integer to binary and pad with zeros
    return [int(x) for x in format(int_representation, '0{}b'.format(num_bits))]

def decode_number(genome:Genome, min_value=-10, max_value=10, num_bits=32):
    # Convert the binary string to an integer
    int_representation = int(''.join(map(str, genome)), 2)
    # Scale down to the normalized value
    normalized = int_representation / (2**num_bits - 1)
    # Denormalize to get the real number
    return min_value + normalized * (max_value - min_value)

def encode_list(number_list,min_value=-10, max_value=10, num_bits=32):
    encoded = []
    for number in number_list:
        encoded += encode_number(number,min_value, max_value, num_bits)
    return encoded

def decode_list(encoded_list,min_value=-10, max_value=10, num_bits=32):
    numbers = []
    for i in range(0, len(encoded_list), num_bits):
        binary_list = encoded_list[i:i + num_bits]
        number = decode_number(binary_list, min_value, max_value, num_bits)
        numbers.append(number)
    return numbers

def generate_population(set_of_points:np.ndarray,min_value=-10, max_value=10, num_bits=32) -> Population:
    return [encode_list(set_of_points[point],min_value, max_value, num_bits) for point in range(len(set_of_points))]

def single_point_crossover(a: Genome, b: Genome, print_cutoff=False) -> Tuple[Genome, Genome]:
    if len(a) != len(b):
        raise ValueError("Genomes a and b must be of same length")
    length = len(a)
    # if the length less than 2, then there is no point to do the function
    if length < 2:  
        return a, b
    # generate random number as the cutoff of the crossover
    p = randint(1, length - 1)
    if print_cutoff == True:
        print(p)
    
    return a[0:p] + b[p:], b[0:p] + a[p:]

def mutation(genome: Genome, num: int = 1, probability: float = 0.5) -> Genome:
    # num: generate how many chromosome(s) that we want to mutate
    for _ in range(num):
        # index sets which chromosome we want to change
        index = randrange(len(genome))
        # the change algorithm
        genome[index] = genome[index] if random() > probability else abs(genome[index] - 1)
    return genome

def fitness_function(genome: Genome, objective_function,min_value=-10, max_value=10, num_bits=32) -> float:
    X = decode_list(genome, min_value, max_value, num_bits)
    return objective_function(X)

# for convenience, call fitness function from these functions below using partial(...)
def population_fitness(population: Population, fitness_func: FitnessFunc) -> float:
    return sum([fitness_func(genome) for genome in population])

def selection_pair(population: Population, fitness_func: FitnessFunc) -> Population:
    return choices(
        population=population,
        weights=[fitness_func(gene) for gene in population],
        k=2
    )
def sort_population(population: Population, fitness_func: FitnessFunc, minimize=True) -> Population:
    return sorted(population,key= lambda x: fitness_func(x), reverse=not minimize)

def genome_to_string(genome: Genome) -> str:
    return "".join(map(str, genome))

def print_stats(population: Population, generation_id: int, fitness_func: FitnessFunc,num_bits=64,binary_mode=False):
    print("GENERATION %02d" % generation_id)
    print("=============")
    if binary_mode == True:
        print("Population: [%s]" % ", ".join([genome_to_string(gene) for gene in population]))
        print("Avg. Fitness: %f" % (population_fitness(population, fitness_func) / len(population)))
        sorted_population = sort_population(population, fitness_func)
        print(
            "Best: %s (%f)" % (genome_to_string(sorted_population[0]), fitness_func(sorted_population[0])))
        print("Worst: %s (%f)" % (genome_to_string(sorted_population[-1]),
                                fitness_func(sorted_population[-1])))
        print("")
    else:
        print(f"Population: {[decode_list(population[i],num_bits=num_bits) for i in range (len(population))]}")
        print("Avg. Fitness: %f" % (population_fitness(population, fitness_func) / len(population)))
        sorted_population = sort_population(population, fitness_func)
        print(f"Best: {(decode_list(population[0],num_bits=num_bits))}")
        print(f"Worst: {(decode_list(population[1],num_bits=num_bits))}")
        print("")

def run_evolution(
        populate_func,
        fitness_func,
        fitness_limit: int,
        minimize = False,
        selection_func: SelectionFunc = selection_pair,
        crossover_func: CrossoverFunc = single_point_crossover,
        mutation_func: MutationFunc = mutation,
        sort_func = sort_population,
        generation_limit: int = 100,
        num_bits:int=64,
        printer: Optional[PrinterFunc] = None) \
        -> Tuple[Population, int]:
    population = populate_func()

    for i in range(generation_limit):
        population = sort_func(population,fitness_func)

        if printer is not None:
            printer(population, i, fitness_func,num_bits)

        if minimize==True:
            cutoff_criteria = fitness_func(population[0])
        else:
            cutoff_criteria = 1-fitness_func(population[0])

        if cutoff_criteria < fitness_limit:
            break

        next_generation = population[0:2]

        for j in range(int(len(population) / 2) - 1):
            parents = selection_func(population)
            offspring_a, offspring_b = crossover_func(parents[0], parents[1])
            offspring_a = mutation_func(offspring_a)
            offspring_b = mutation_func(offspring_b)
            next_generation += [offspring_a, offspring_b]

        population = next_generation

    return population, i

n_point = 100
k_max = 100
dim = 3
epsilon = 10**(-7)
def objective_function(x):
    f=0
    for i in range (len(x)):
        f += (x[i]**2-10*np.cos(2*np.pi*x[i])+10)
    return f
boundaries = np.array([(-5,5) for _ in range (dim)])
min_value = boundaries.min()
max_value = boundaries.max()
num_bits = 64  # Number of bits for each number
iter_points = generate_points(dim,n_point,boundaries[:,0],boundaries[:,1])
popA = generate_population(iter_points,num_bits=num_bits,min_value=min_value,max_value=max_value)
decoded_genome = [decode_list(popA[i],num_bits=num_bits) for i in range (len(popA))]
population,generation = run_evolution(
    populate_func=partial(
        generate_population,set_of_points = iter_points,num_bits=num_bits,min_value=min_value,max_value=max_value
    ),
    fitness_func= partial(
        fitness_function, objective_function=objective_function,num_bits=num_bits,min_value=min_value,max_value=max_value
    ),
    minimize=True,
    sort_func=partial(
        sort_population,minimize=True),
    selection_func=partial(selection_pair,
                           fitness_func=partial(
                               fitness_function,
                               num_bits=num_bits,
                               min_value=min_value,max_value=max_value
                           )),
    fitness_limit=epsilon,
    generation_limit=k_max
)
print(decode_list(population[0],num_bits=num_bits))