import numpy as np
from functools import partial
from scipy.stats import qmc

class RADE:
    def __init__(
            self,
            system_equation,
            boundaries,
            population_size,
            max_generation,
            number_per_subpopulation,
            theta,
            tau_d,
            F_init,
            CR_init,
            max_archive_size,
            max_memories_size,
            beta = 100,
            rho = 1e-8,
            seed = None
            ) -> None:
        """
        PROPERTIES
            archive: list= list of eligible roots
            systemeq: func= systen of equation in question
            boundaries: np.ndarray= boundaries/constraint of the problem
            dim: int= dimension of roots
            population_size: int= population size per generation
            max_generation: int= maximal iteration
            max_memories_size: int=  maximum memories size
            memories_F: np.ndarray= memories of scaling factor
            memories_CR: np.ndarray= memories of crossover rate
            number_per_subpopulation: int= number of individuals in a subpopulation of the target individu
            theta: the error of the objective function
            tau_d: the minimum distance between roots
            max_archive_size: maximum found roots
            beta: int= large constant to control the penalty scale
            rho: float= small constant to adjust the radius of the repulsion regions
            seed: int= random state
        """
        self.boundaries = boundaries
        self.systemeq = system_equation
        self.population_size = population_size
        self.max_generation = max_generation
        self.number_per_subpopulation = number_per_subpopulation
        self.tau_d = tau_d
        self.max_archive_size = max_archive_size
        self.theta = theta
        self.seed = seed
        self.dim = boundaries.shape[0]
        self.archive = []
        self.max_memories_size = max_memories_size
        self.memories_F = np.ones(max_memories_size)*F_init
        self.memories_CR = np.ones(max_memories_size)*CR_init
        self.beta = beta
        self.rho = rho
        return None
    
    def objective_function(self,x:np.ndarray):
        res = 0
        F_array = self.systemeq(x)
        for f in F_array:
            res += np.abs(f)
        self.condition = -1+self.epsilon
        return -1/(1+res)
    
    def chi_p(self, delta_j, rho):
        """
        Characteristic Function
        """
        return 1 if delta_j <= rho else 0
    
    def repulsion_function(
            self,
            x):
        """
        Repulsion Function
        """
        f_x = self.objective_function(x)
        Rx = 0
        for x_star in self.archive:
            delta_j = np.linalg.norm(x-x_star)
            Rx += np.exp(-delta_j) * self.chi_p(delta_j, self.rho)
        Rx *= self.beta
        Rx += f_x
        return Rx
    
    def fitness_function(
            self,
            x
            ):
        """
        Fitness Function
        """ 
        if self.archive == []:
            return self.objective_function(x)
        else:
            return self.repulsion_function(x)
        
        """GENERATE POINTS USING SOBOL SEQUENCE"""
    def generate_points(
            self,
            npoint: int, 
            low=-10, 
            high=10, 
            sobol=True
            ):
        """
                Generates points within the specified bounds.

            Args:
                dim: Number of dimensions.
                npoint: Number of points to generate.
                low: Lower bound for each variable (scalar or list/numpy array).
                high: Upper bound for each variable (scalar or list/numpy array).
                sobol: Flag indicating whether to use Sobol sequence (True) or random sampling (False).

            Returns:
                A numpy array of size (npoint, dim) representing the generated points.
        """

        if type(low) != type(high):
            raise TypeError('The type of "low" and "high" should be the same.')

        # Handle boundaries
        if type(low) == int:
            boundaries = [(low, high) for _ in range(self.dim)]
        elif type(low) in (list, np.ndarray):
            if len(low) != len(high):
                raise TypeError('The length of "low" and "high" should be the same.')
            else:
                boundaries = [(low[i], high[i]) for i in range(len(low))]

        # Generate points based on the sobol flag
        if sobol:
            sampler = qmc.Sobol(d=self.dim,scramble=True,seed=self.seed)
            sample = sampler.random(n=npoint)
            scaled_points = qmc.scale(sample=sample,l_bounds=low,u_bounds=high)

        else:
            # Generate random points
            np.random.seed(self.seed)
            scaled_points = np.zeros((npoint, self.dim))
            for i in range(self.dim):
                min_val, max_val = boundaries[i]
                scaled_points[:, i] = np.random.uniform(min_val, max_val, npoint)

        return scaled_points
    
    def mutate(self,population, F):
        """
        Mutation function for DE
        """
        # Vectorized mutation operation
        indices = np.arange(population.shape[0])
        np.random.shuffle(indices)
        r1, r2, r3 = population[indices[:3]]
        mutant = r1 + F * (r2 - r3)
        return mutant
    
    def crossover(self, target, mutant, CR):
        cross_points = np.random.random(target.shape) < CR
        # Ensure at least one true crossover point
        if not np.any(cross_points):
            cross_points[np.random.randint(0, target.shape[0])] = True
        trial = np.where(cross_points, mutant, target)
        return trial
    
    def mutation_penalty(
            self,
            x_i, 
            subpop_i, 
            boundaries, 
            scaling_factor
            ):
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
        dv_i = self.mutate(subpop_i_copy,scaling_factor)

        """Set penalty for every donor vector that violates the boundaries"""
        for j in range (len(dv_i)):
            if dv_i[j] < boundaries[j,0]:
                dv_i[j] = (x_i[j]+boundaries[j,0])/2
            elif dv_i[j] > boundaries[j,1]:
                dv_i[j] = (x_i[j]+boundaries[j,1])/2
        return dv_i
    
    def subpopulating(
            self,
            individual, 
            population, 
            t,
            return_index = False,
            show_distances = False
            ): 
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
            
    def update_archive(
            self,
            x: np.ndarray
            ):
        """
        INPUT:
            x: list= Individual

        OUTPUT: 
            archive: list= updated archive
        """

        f_x = self.objective_function(x)
        s = len(self.archive) # archive current size
        if f_x < self.theta: # x is a root
            # print(f'f({x})= {f_x}')
            if s == 0: # archive is empty
                self.archive.append(x)
                s+=1
            else:
                """Find the closest solution x_prime ∈ archive to x in the decision space"""
                dist_min = np.linalg.norm(x-self.archive[0])
                idx_min = 0
                x_prime= self.archive[idx_min]
                for i in range(1,len(self.archive)): 
                    dist = np.linalg.norm(x-self.archive[i])
                    if dist < dist_min:
                        dist_min = dist
                        x_prime = self.archive[i]
                        idx_min = i
                f_x_prime = self.objective_function(x_prime)
                if dist_min < self.tau_d: # x and x_prime are too close
                    if f_x < f_x_prime:
                        x_prime = x
                        self.archive[idx_min] = x_prime
                else:
                    if s < self.max_archive_size:
                        self.archive.append(x)
                        s += 1
                    else:       # archive is full
                        if f_x<f_x_prime:
                            x_prime = x
                            self.archive[idx_min] = x_prime
    
    def update_parameter(self):
        """
        OUTPUT
            Fi: float= scaling factor
            CRi: float= crossover rate
        """
        # Randomly select an index
        hi = np.random.randint(0, self.max_memories_size)
        # Generate Fi using the Cauchy distribution with the location parameter MF[hi] and scale 0.1
        Fi = np.random.standard_cauchy() * 0.1 + self.memories_F[hi]
        # Generate CRi using the Gaussian distribution with mean MCR[hi] and standard deviation 0.1
        CRi = np.random.normal(self.memories_CR[hi], 0.1)
        # Ensure CRi is within the range [0, 1] and Fi is within the range [0,2]
        Fi = np.clip(Fi, 0, 1)
        CRi = np.clip(CRi, 0, 1)
        return Fi, CRi
    
    def meanWL(self, elements, weights):
        """
        Calculate the weighted Lehmer mean of elements.
        Lehmer mean is calculated as the weighted sum of the squares
        divided by the weighted sum of the elements.
        """
        numerator = np.sum(np.multiply(np.square(elements), weights))
        denominator = np.sum(np.multiply(elements, weights))
        return numerator / denominator if denominator != 0 else 0
    
    # Define the weighted arithmetic mean function
    def meanWA(self, elements, weights):
        """
        Calculate the weighted arithmetic mean of elements.
        This is the standard weighted mean.
        """
        return np.average(elements, weights=weights)
    
    def update_history(
            self,
            S_F,
            S_CR,
            k
            ):
        weights = np.array([1 for _ in range (len(S_F))])
        if len(S_F)!=0:
            self.memories_F[k] = self.meanWL(S_F,weights) 
        if len(S_CR)!=0:
            self.memories_CR[k] = self.meanWA(S_CR,weights)
    
    def RADE_evaluation(
            self,
            print_gen = False,
            root_history = False
            ):
        """
        OUTPUT
            best: np.ndarray= best individual in term of objective function
            fitness[best_idx]: np.ndarray= best value of objective function
            archive: roots found in the problem
        """
        np.random.seed(self.seed)
        population = self.generate_points(
            npoint=self.population_size,
            low=self.boundaries[:,0],
            high=self.boundaries[:,1],
            sobol=True)
        fitness = np.asarray([self.objective_function(ind) for ind in population])
        best_idx = np.argmin(fitness)
        best = population[best_idx] # objective function must be minimizing
        subpop = np.array([self.subpopulating(individual=xi, population=population, t=self.number_per_subpopulation) for xi in population])

        k=0
        if root_history == True:
            history_root = {}
        for gen in range(self.max_generation):
            S_F, S_CR = [],[]
            for i in range(self.population_size):
                F_i,CR_i = self.update_parameter()
                x_i = population[i]
                mutant = self.mutation_penalty(
                    x_i=x_i,
                    subpop_i=subpop[i],
                    boundaries=self.boundaries,
                    scaling_factor=F_i)
                trial = self.crossover(target=population[i], mutant=mutant, CR=CR_i)
                trial_fitness = self.fitness_function(x=trial)
                
                if trial_fitness < fitness[i]:
                    fitness[i] = trial_fitness
                    population[i] = trial
                    self.update_archive(x=trial)
                    S_F.append(F_i)
                    S_CR.append(CR_i)
                    if trial_fitness < fitness[best_idx]:
                        best_idx = i
                        best = trial

            if print_gen == True:
                print(f"=========Generation {gen}=========")
                # print(f"Best Fitness: {fitness[best_idx]}")
                print(f'Archive:{self.archive}')
            
            if root_history == True:
                history_root[gen] = self.archive.copy()
            
                # print(f'S_F: {S_F}\nS_CR: {S_CR}')
            """Update parameter history"""
            if (len(S_F)!=0) & (len(S_CR)!=0):
                self.update_history(S_F=S_F,S_CR=S_CR,k=k)
                # print(f'M_F: {M_F}\nM_CR: {M_CR}')
                k +=1
                if k >= self.max_memories_size:
                    k = 1

        
        if root_history == True:
            return best, fitness[best_idx], self.archive, history_root
        else:
            return best, fitness[best_idx], self.archive