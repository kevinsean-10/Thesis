import numpy as np
from functools import partial
import sobol_seq

"""GENERATE POINTS USING SOBOL SEQUENCE"""
def generate_points(dim: int,
                    npoint:int,
                    low=-10,
                    high=10,
                    sobol = True):
    if type(low) != type(high):
        raise TypeError('The type of "low" and "high" should be the same.')
    if type(low) == int:
        boundaries = [(low,high) for _ in range (dim)]
    elif type(low) == list or type(low) == np.ndarray:
        if len(low) != len(high):
            raise TypeError('The length of "low" and "high" should be the same.')
        else:
            boundaries = [(low[i],high[i]) for i in range (len(low))]

    if sobol == True:
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
    
    else:
        scaled_points = np.zeros((npoint, dim))
        for i in range(dim):
            min_val, max_val = boundaries[i]
            scaled_points[:, i] = np.random.uniform(min_val, max_val, npoint)

    return scaled_points

"""Roullete wheel selection"""
def selection(population: np.ndarray,
              fitness: np.ndarray):
    population_size = population.shape[0]
    selection_probs = np.array([1 / (fit + 1) for fit in fitness]) # add one to avoid negative probability
    total_probs = sum(selection_probs)
    selection_probs = np.array([prob / total_probs for prob in selection_probs])
    selected_indices = np.random.choice(a=np.arange(population_size),size=population_size,p=selection_probs)
    selected_population = np.array([population[i] for i in selected_indices])
    return selected_population

"""Single point crossover"""
def crossover(parent1, parent2):
    dimension = len(parent1)
    crossover_point = np.random.randint(1, dimension)
    offspring1 = np.append(parent1[:crossover_point], parent2[crossover_point:])
    offspring2 = np.append(parent2[:crossover_point], parent1[crossover_point:])
    return [offspring1, offspring2]

"""One point mutation"""
def mutate(individual,mutation_rate, boundaries):
    for j in range(len(individual)):
        if np.random.random() < mutation_rate:
            individual[j] = np.random.uniform(boundaries[j][0], boundaries[j][1])
    return individual


def recombination(population: np.ndarray,mutation_rate, boundaries):
    offspring_population = []
    population_size = population.shape[0]
    for i in range(0, population_size, 2):
        parent1 = population[i]
        parent2 = population[i + 1]
        offspring1, offspring2 = crossover(parent1=parent1,parent2=parent2)
        offspring1 = mutate(individual=offspring1,mutation_rate=mutation_rate, boundaries=boundaries)
        offspring2 = mutate(individual=offspring2,mutation_rate=mutation_rate, boundaries=boundaries)
        offspring_population.extend([offspring1, offspring2])
    offspring_population = np.array(offspring_population)
    return offspring_population

def GA(objective_function,
       population_size,
       boundaries,
       max_generation,
       mutation_rate,
       seed=0,
       print_stat = False):
    
    np.random.seed(seed)
    dim = boundaries.shape[1]
    population = generate_points(dim=dim,npoint=population_size,low=boundaries[:,0],high=boundaries[:,1],sobol=True)
    fitness = np.asarray([objective_function(ind) for ind in population])
    best_idx = np.argmin(fitness)
    best_individual = population[best_idx]

    for generation in range(max_generation):
        selected_population = selection(population=population,fitness=fitness)
        offspring_population = recombination(population=selected_population,
                                             mutation_rate=mutation_rate,
                                             boundaries=boundaries)
        population = offspring_population
        fitness = np.asarray([objective_function(ind) for ind in population])
        best_idx = np.argmin(fitness)
        best_individual = population[best_idx]
        best_fitness = fitness[best_idx]
        if print_stat == True:

            print(f"=========Generation {generation}=========")
            print(f"Best Individual: {best_individual}")
            print(f"Best Score: {best_fitness}\n")


    return best_individual, best_fitness

