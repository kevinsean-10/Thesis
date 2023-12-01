import numpy as np
import sobol_seq
import matplotlib.pyplot as plt
from functools import partial

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


def crossover(individual:np.ndarray, 
              honor_vector:np.ndarray,
              crossover_rate,
              dimension):
    trial_vector = np.zeros(shape=individual.shape)
    for j in range (individual.shape[0]):
        rand_j = np.random.random() # a uniformly distributed random number from [0, 1]
        j_rand = np.random.randint(dimension+1) #a random integer uniformly generated from {1, . . . , n}
        if (rand_j<crossover_rate) | (j == j_rand):
            trial_vector[j] = honor_vector[j]
        else: #(rand_j>=CR) and (j != j_rand)
            trial_vector[j] = individual[j]
    return trial_vector

def mutation(xr1,xr2,xr3,mutation_factor):
    donor_vector = xr1 + mutation_factor*(xr2-xr3)
    return donor_vector

def reproduction(population,
                 objective_function,
                 dimension,
                 mutation_factor,
                 crossover_rate,
                 seed=0):
    np.random.seed(seed=seed)
    pop_ids = np.arange(population.shape[0])
    for i in range (len(population)):
        x_i = population[i]
        """GENERATE: three distinct individuals xr1, xr1, xr1 from the current population randomly"""
        pop_ids_no_i = np.delete(pop_ids,i)
        xr1,xr2,xr3 = population[np.random.choice(pop_ids_no_i,3,replace=False)] 

        """MUTATION: Form the donor/mutation vector"""
        dv_i = mutation(xr1,xr2,xr3,mutation_factor)

        """CROSSOVER: The trial vector ui is developed either from the elements of the target vector xi or the elements of the donor vector vi"""
        tv_i = crossover(x_i,dv_i,crossover_rate,dimension)

        """EVALUATE: If f (ui) ≤ f (xi) then replace the individual xi in the population with the trial vector ui"""
        if objective_function(tv_i)<=objective_function(x_i):
            population[i] = tv_i
    population = np.array(sorted(population,key=lambda x: objective_function(x)))
    return population

def differensial_evolution(number_of_population,
                           objective_function,
                           boundaries, 
                           dimension, 
                           gen_max, 
                           mutation_factor, 
                           crossover_rate,
                           seed=0):

    popA = generate_points(dim=dimension,
                            npoint=number_of_population,
                            low=boundaries[:,0],
                            high=boundaries[:,1])
    genA = {0:popA}
    # print(genA)
    generation = 0
    while generation<gen_max:
        popB = genA[generation].copy()
        popB = reproduction(population=popB,
                            objective_function=objective_function,
                            dimension=dimension,
                            mutation_factor=mutation_factor,
                            crossover_rate=crossover_rate,
                            seed=seed)
        generation += 1
        genA[generation] = popB
    return genA
    # print(genA[FEs_max/NP])
