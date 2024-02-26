import numpy as np
import sobol_seq
import matplotlib.pyplot as plt
from functools import partial

def generate_points(dim: int,
                    npoint:int,
                    low=-10,
                    high=10,
                    sobol = True):
    """
    Generate points with option to use sobol sequence
    """
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

def update_archive(x: np.ndarray,
                   objective_function,
                   archive,
                   theta,
                   tau_d,
                   s_max):
    """
    INPUT:
        x: list= Individual
        theta: float/int= accuracy level
        tau_d: float= distance radius
        s_max: int: maximum archive size
        archive: list= archive
        s : int= archive current size

    OUTPUT: 
        archive: list= updated archive
    """

    f_x = objective_function(x)
    s = len(archive) # archive current size
    if f_x < theta: # x is a root
        # print(f'f({x})= {f_x}')
        if s == 0: # archive is empty
            archive.append(x)
            s+=1
        else:
            """Find the closest solution x_prime ∈ archive to x in the decision space"""
            dist_min = np.linalg.norm(x-archive[0])
            idx_min = 0
            x_prime= archive[idx_min]
            for i in range(1,len(archive)): 
                dist = np.linalg.norm(x-archive[i])
                if dist < dist_min:
                    dist_min = dist
                    x_prime = archive[i]
                    idx_min = i
            f_x_prime = objective_function(x_prime)
            if dist_min < tau_d: # x and x_prime are too close
                if f_x < f_x_prime:
                    x_prime = x
                    archive[idx_min] = x_prime
            else:
                if s < s_max:
                    archive.append(x)
                    s += 1
                else:       # archive is full
                    if f_x<f_x_prime:
                        x_prime = x
                        archive[idx_min] = x_prime
    return archive

def crossover(individual:np.ndarray, 
              honor_vector:np.ndarray,
              crossover_rate):
    trial_vector = np.zeros(shape=individual.shape)
    for j in range (individual.shape[0]):
        rand_j = np.random.random() # a uniformly distributed random number from [0, 1]
        j_rand = np.random.randint(individual.shape[0]+1) #a random integer uniformly generated from {1, . . . , n}
        if (rand_j<crossover_rate) | (j == j_rand):
            trial_vector[j] = honor_vector[j]
        else: #(rand_j>=CR) and (j != j_rand)
            trial_vector[j] = individual[j]
    return trial_vector

def mutate(population, mutation_factor):
    """
    Mutation function for DE
    """
    # Vectorized mutation operation
    indices = np.arange(population.shape[0])
    np.random.shuffle(indices)
    r1, r2, r3 = population[indices[:3]]
    mutant = r1 + mutation_factor * (r2 - r3)
    return mutant

def mutation_penalty(x_i, 
                     population, 
                     boundaries, 
                     mutation_factor):
    """
    INPUT:
        x_i: np.ndarray= target x_i
        population: np.ndarray= number of individuals closest to x_i
        boundaries: np.ndarray= boundaries/constraints of the function
        scaling_factor: scaling factor of the function

    OUTPUT:
        dv_i: np.ndarray= donor vector that has been mutated and penalized.
    """

    """Generate three distinct individuals xr1, xr1, xr1 from the current population randomly"""
    population_copy = population.copy()
    pop_ids = np.arange(population_copy.shape[0])
    indices_to_delete = np.where(np.all(population_copy == x_i, axis=1))[0] # Ensure that x_i is excluded from the selected subpopulation
    pop_ids_no_i = np.delete(pop_ids, indices_to_delete, axis=0)
    population_copy = population_copy[pop_ids_no_i]

    """Mutation form the donor/mutation vector"""
    dv_i = mutate(population_copy,mutation_factor)

    """Set penalty for every donor vector that violates the boundaries"""
    for j in range (len(dv_i)):
        if dv_i[j] < boundaries[j,0]:
            dv_i[j] = (x_i[j]+boundaries[j,0])/2
        elif dv_i[j] > boundaries[j,1]:
            dv_i[j] = (x_i[j]+boundaries[j,1])/2
    return dv_i

def reproduction(population,
                 objective_func,
                 boundaries,
                 mutation_factor,
                 crossover_rate,
                 seed=None):
    np.random.seed(seed=seed)
    for i in range (len(population)):
        x_i = population[i]

        """MUTATION: Form the donor/mutation vector"""
        dv_i = mutation_penalty(x_i=x_i,
                                population=population,
                                boundaries=boundaries,
                                mutation_factor=mutation_factor)

        """CROSSOVER: The trial vector ui is developed either from the elements of the target vector xi or the elements of the donor vector vi"""
        tv_i = crossover(x_i,dv_i,crossover_rate)

        """EVALUATE: If f (ui) ≤ f (xi) then replace the individual xi in the population with the trial vector ui"""
        if objective_func(tv_i)<=objective_func(x_i):
            population[i] = tv_i

    # Sort the population based on the lowest score (best) on the objective function
    population = np.array(sorted(population,key=lambda x: objective_func(x)))
    return population

# def differensial_evolution(number_of_population,
#                            objective_func,
#                            population_size,
#                            boundaries, 
#                            dimension, 
#                            max_generation, 
#                            mutation_factor, 
#                            crossover_rate,
#                            seed=0,
#                            history = False):
#     np.random.seed(seed)
#     population = generate_points(dim=dimension,
#                                  npoint=population_size,
#                                  low=boundaries[:,0],
#                                  high=boundaries[:,1],
#                                  sobol=True)
#     fitness = np.asarray([objective_func(ind) for ind in population])
#     best_idx = np.argmin(fitness)
#     best = population[best_idx] # objective function must be minimizing

#     if history == True:
#         hist = {0:population}

#     for gen in range (max_generation):
#         popB = genA[generation].copy()
#         popB = reproduction(population=popB,
#                             objective_function=objective_function,
#                             dimension=dimension,
#                             mutation_factor=mutation_factor,
#                             crossover_rate=crossover_rate,
#                             seed=seed)
#         if history == True:
#             hist[gen+1] = popB
#     return genA
#     # print(genA[FEs_max/NP])

def differensial_evolution(objective_func,
                           population_size,
                           boundaries, 
                           gen_max, 
                           mutation_factor, 
                           crossover_rate,
                           seed=None,
                           print_gen = False,
                           history = False):
    np.random.seed(seed)
    dimension = len(boundaries)
    population = generate_points(dim=dimension,
                                 npoint=population_size,
                                 low=boundaries[:,0],
                                 high=boundaries[:,1],
                                 sobol=True)
    fitness = np.asarray([objective_func(ind) for ind in population])
    best_idx = np.argmin(fitness)
    best = population[best_idx] # objective function must be minimizing

    if history == True:
        hist = {0:population}

    for gen in range (gen_max):

        for i in range (population_size):
            x_i = population[i]
            dv_i = mutation_penalty(x_i=x_i,
                                    population=population,
                                    boundaries=boundaries,
                                    mutation_factor=mutation_factor)
            trial = crossover(x_i,dv_i,crossover_rate)
            trial_fitness = objective_func(trial)

            if trial_fitness<=fitness[i]:
                fitness[i] = objective_func(trial)
                population[i] = trial                
                if trial_fitness < fitness[best_idx]:
                    best_idx = i
                    best = trial

        if print_gen == True:
            print(f"=========Generation {gen}=========")
            # print(f"Best Fitness: {fitness[best_idx]}")
            print(f'Best Point:{best} with score {fitness[best_idx]}')       

        if history == True:
            hist[gen+1] = population
            print(hist)

    if history == True:
        return best, fitness[best_idx]
    else:
        return best, fitness[best_idx]
