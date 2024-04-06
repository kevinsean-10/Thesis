import numpy as np
from functools import partial
import sobol_seq

def chi_p(delta_j, rho):
    """
    Characteristic Function
    """
    return 1 if delta_j <= rho else 0

def repulsion_function(x,
                       archive,
                       objective_func,
                       beta=100,
                       rho=1e-8):
    """
    Repulsion Function
    """
    f_x = objective_func(x)
    Rx = 0
    for x_star in archive:
        delta_j = np.linalg.norm(x-x_star)
        Rx += np.exp(-delta_j) * chi_p(delta_j, rho)
    Rx *= beta
    Rx += f_x
    return Rx

def fitness_function(x,
                     archive,
                     objective_func,
                     repulsion_func=repulsion_function):
    """
    Fitness Function
    """
    f_x = objective_func(x)
    if archive == []:
        return f_x
    else:
        return repulsion_func(x,archive)
    
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



def mutate(population, F):
    """
    Mutation function for DE
    """
    # Vectorized mutation operation
    indices = np.arange(population.shape[0])
    np.random.shuffle(indices)
    r1, r2, r3 = population[indices[:3]]
    mutant = r1 + F * (r2 - r3)
    return mutant

def crossover(target, mutant, CR):
    cross_points = np.random.random(target.shape) < CR
    # Ensure at least one true crossover point
    if not np.any(cross_points):
        cross_points[np.random.randint(0, target.shape[0])] = True
    trial = np.where(cross_points, mutant, target)
    return trial

def mutation_penalty(x_i, 
                     subpop_i, 
                     boundaries, 
                     scaling_factor):
    """
    INPUT:
        x_i: np.ndarray= target x_i
        subpop_i: np.ndarray= number of individuals closest to x_i
        boundaries: np.ndarray= boundaries/constraints of the function
        scaling_factor: scaling factor of the function

    OUTPUT:
        dv_i: np.ndarray= donor vector that has been mutated and penalized.
    """

    """Generate three distinct individuals xr1, xr1, xr1 from the current population randomly"""
    subpop_i_copy = subpop_i.copy()
    pop_ids = np.arange(subpop_i_copy.shape[0])
    indices_to_delete = np.where(np.all(subpop_i_copy == x_i, axis=1))[0] # Ensure that x_i is excluded from the selected subpopulation
    subpop_ids_no_i = np.delete(pop_ids, indices_to_delete, axis=0)
    subpop_i_copy = subpop_i_copy[subpop_ids_no_i]

    """Mutation form the donor/mutation vector"""
    dv_i = mutate(subpop_i_copy,scaling_factor)

    """Set penalty for every donor vector that violates the boundaries"""
    for j in range (len(dv_i)):
        if dv_i[j] < boundaries[j,0]:
            dv_i[j] = (x_i[j]+boundaries[j,0])/2
        elif dv_i[j] > boundaries[j,1]:
            dv_i[j] = (x_i[j]+boundaries[j,1])/2
    return dv_i

def subpopulating(individual, 
                  population, 
                  t,
                  return_index = False,
                  show_distances = False): 
    """
    Calculate Euclidean distances and select t closest individuals.
    INPUT:
        individual: np.ndarray= individual target
        population: np.ndarray= population where the target in
        t: int= max number of units in a subpopulation

    OUTPUT:
        closest_indices: np.ndarray= array of closest index of the target individual
        subpop: np.ndarray= array of closest individuals of the target individual
    """
    
    """Calculate the Euclidean distances from the individual to all others in the population"""
    distances = np.sqrt(np.sum((population - individual) ** 2, axis=1))

    # Get the indices of the individuals with the smallest distances
    closest_indices = np.argsort(distances)[:t]
    # Form the subpopulation with the closest individuals
    subpop = population[closest_indices]

    if show_distances == True:
        print(f'Distance: \n{distances[:t]}')
    if return_index == True:
        if t == 1:
            return closest_indices,subpop.flatten()
        else:
            return closest_indices,subpop
    else:
        if t == 1:
            return subpop.flatten()
        else:
            return subpop

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

def update_parameter(M_F,
                     M_CR,
                     Hm:int):
    """
    INPUT
        M_F: np.ndarray= Historical memories of scaling factor of DE as F
        M_CR: np.ndarray= Historical memories crossover rate of DE as CR
        Hm: int= Size of Historical Memories
    
    OUTPUT
        Fi: float= scaling factor
        CRi: float= crossover rate
    """
    # Randomly select an index
    hi = np.random.randint(0, Hm)
    # Generate Fi using the Cauchy distribution with the location parameter MF[hi] and scale 0.1
    Fi = np.random.standard_cauchy() * 0.1 + M_F[hi]
    # Generate CRi using the Gaussian distribution with mean MCR[hi] and standard deviation 0.1
    CRi = np.random.normal(M_CR[hi], 0.1)
    # Ensure CRi is within the range [0, 1] and Fi is within the range [0,2]
    Fi = np.clip(Fi, 0, 1)
    CRi = np.clip(CRi, 0, 1)
    return Fi, CRi

def meanWL(elements, weights):
    """
    Calculate the weighted Lehmer mean of elements.
    Lehmer mean is calculated as the weighted sum of the squares
    divided by the weighted sum of the elements.
    """
    numerator = np.sum(np.multiply(np.square(elements), weights))
    denominator = np.sum(np.multiply(elements, weights))
    return numerator / denominator if denominator != 0 else 0

# Define the weighted arithmetic mean function
def meanWA(elements, weights):
    """
    Calculate the weighted arithmetic mean of elements.
    This is the standard weighted mean.
    """
    return np.average(elements, weights=weights)

def update_history(M_F,M_CR,S_F,S_CR,k):
    weights = np.array([1 for _ in range (len(S_F))])
    if len(S_F)!=0:
        M_F[k] = meanWL(S_F,weights) 
    if len(S_CR)!=0:
        M_CR[k] = meanWA(S_CR,weights)
    return M_F,M_CR

def RADE(archive, 
         objective_func,
         bounds, 
         population_size, 
         max_generation, 
         memories_F, 
         memories_CR,
         num_l,
         theta,
         tau_d,
         archive_size_max,
         seed=0,
         print_gen = False,
         root_history = False):
    """
    INPUT
        archive: list= list of eligible roots
        objective_func: func= objective function to solve the NLS
        bounds: np.ndarray= boundaries/constraint of the problem
        population_size: int= population size per generation
        max_generation: int= maximal iteration
        memories_F: np.ndarray= memories of scaling factor
        memories_CR: np.ndarray= memories of crossover rate
        num_l: int= number of individuals in a subpopulation of the target individu
        theta: the error of the objective function
        tau_d: the minimum distance between roots
        archive_size_max: maximum found roots

    OUTPUT
        best: np.ndarray= best individual in term of objective function
        fitness[best_idx]: np.ndarray= best value of objective function
        archive: roots found in the problem
    """
    np.random.seed(seed)
    dimensions = len(bounds)
    population = generate_points(dim=dimensions,
                                 npoint=population_size,
                                 low=bounds[:,0],
                                 high=bounds[:,1],
                                 sobol=True)
    fitness = np.asarray([objective_func(ind) for ind in population])
    best_idx = np.argmin(fitness)
    best = population[best_idx] # objective function must be minimizing
    subpopA = np.array([subpopulating(xi, population, num_l) for xi in population])
    Hm = len(memories_F)
    k=0

    if root_history == True:
        history_root = {}
    for gen in range(max_generation):
        S_F, S_CR = [],[]
        for i in range(population_size):
            F_i,CR_i = update_parameter(memories_F,memories_CR,Hm)
            x_i = population[i]
            mutant = mutation_penalty(x_i,subpopA[i],bounds,F_i)
            trial = crossover(population[i], mutant, CR_i)
            trial_fitness = fitness_function(trial, 
                                             archive,
                                             objective_func=objective_func,
                                             repulsion_func=partial(repulsion_function,
                                                                    objective_func = objective_func))
            
            if trial_fitness < fitness[i]:
                fitness[i] = trial_fitness
                population[i] = trial
                archive = update_archive(x = trial,
                                        objective_function=objective_func,
                                        archive=archive,
                                        theta=theta,
                                        tau_d=tau_d,
                                        s_max=archive_size_max)
                S_F.append(F_i)
                S_CR.append(CR_i)
                if trial_fitness < fitness[best_idx]:
                    best_idx = i
                    best = trial

        if print_gen == True:
            print(f"=========Generation {gen}=========")
            # print(f"Best Fitness: {fitness[best_idx]}")
            print(f'Archive:{archive}')
        
        if root_history == True:
            history_root[gen] = archive.copy()
        
            # print(f'S_F: {S_F}\nS_CR: {S_CR}')
        """Update parameter history"""
        if (len(S_F)!=0) & (len(S_CR)!=0):
            memories_F,memories_CR = update_history(memories_F,memories_CR,S_F,S_CR,k)
            # print(f'M_F: {M_F}\nM_CR: {M_CR}')
            k +=1
            if k >= Hm:
                k = 1

    
    if root_history == True:
        return best, fitness[best_idx], archive, history_root
    else:
        return best, fitness[best_idx], archive