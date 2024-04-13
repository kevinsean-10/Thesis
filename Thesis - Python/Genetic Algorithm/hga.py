from itertools import product,combinations
import numpy as np
from scipy.stats import qmc

class HGA:
    def __init__(
            self,
            systemeq,
            population_size,
            boundaries,
            max_generation,
            mutation_rate,
            parts,
            epsilon,
            delta,
            seed
            ) -> None:
        """
        PROPERTIES
            systemeq: func= systen of equation in question
            boundaries: np.ndarray= boundaries/constraint of the problem
            dim: int= dimension of roots
            population_size: int= population size per generation
            max_generation: int= maximal iteration
            max_memories_size: int=  maximum memories size
            mutation_rate: np.ndarray= probability of mutation
            epsilon: the error of the objective function
            delta: the minimum distance between roots
            seed: int= random state
        """
        self.population_size = population_size
        self.systemeq = systemeq
        self.boundaries = boundaries
        self.dim = boundaries.shape[0]
        self.max_generation = max_generation
        self.mutation_rate = mutation_rate
        self.parts = parts
        self.epsilon = epsilon
        self.delta = delta
        self.seed = seed
        return None
    
    def objective_function(self,x:np.ndarray):
        res = 0
        F_array = self.systemeq(x)
        for f in F_array:
            res += np.abs(f)
        self.condition = -1+self.epsilon
        return -1/(1+res)
    
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
    
    """Roullete wheel selection"""
    def selection(
            self,
            population: np.ndarray,
            fitness: np.ndarray):
        population_size = population.shape[0]
        selection_probs = np.array([1 / (fit + 1) for fit in fitness]) # add one to avoid negative probability
        total_probs = sum(selection_probs)
        selection_probs = np.array([prob / total_probs for prob in selection_probs])
        selected_indices = np.random.choice(a=np.arange(population_size),size=population_size,p=selection_probs)
        selected_population = np.array([population[i] for i in selected_indices])
        return selected_population
    
    """Single point crossover"""
    def crossover(self,parent1, parent2):
        dimension = len(parent1)
        crossover_point = np.random.randint(1, dimension)
        offspring1 = np.append(parent1[:crossover_point], parent2[crossover_point:])
        offspring2 = np.append(parent2[:crossover_point], parent1[crossover_point:])
        return [offspring1, offspring2]
    
    """One point mutation"""
    def mutate(self, individual, mutation_rate, boundaries):
        for j in range(len(individual)):
            if np.random.random() < mutation_rate:
                individual[j] = np.random.uniform(boundaries[j][0], boundaries[j][1])
        return individual
    
    def recombination(self,population: np.ndarray,mutation_rate, boundaries):
        offspring_population = []
        population_size = population.shape[0]
        for i in range(0, population_size, 2):
            parent1 = population[i]
            parent2 = population[i + 1]
            offspring1, offspring2 = self.crossover(parent1=parent1,parent2=parent2)
            offspring1 = self.mutate(individual=offspring1,mutation_rate=mutation_rate, boundaries=boundaries)
            offspring2 = self.mutate(individual=offspring2,mutation_rate=mutation_rate, boundaries=boundaries)
            offspring_population.extend([offspring1, offspring2])
        offspring_population = np.array(offspring_population)
        return offspring_population
    
    def GA(self,
           population_size,
           boundaries,
           max_generation,
           mutation_rate,
           seed=0,
           print_stat = False):
        
        np.random.seed(seed)
        dim = boundaries.shape[1]
        population = self.generate_points(npoint=population_size,low=boundaries[:,0],high=boundaries[:,1],sobol=True)
        fitness = np.asarray([self.objective_function(ind) for ind in population])
        best_idx = np.argmin(fitness)
        best_individual = population[best_idx]

        for generation in range(max_generation):
            selected_population = self.selection(population=population,fitness=fitness)
            offspring_population = self.recombination(population=selected_population,
                                                      mutation_rate=mutation_rate,
                                                      boundaries=boundaries)
            population = offspring_population
            fitness = np.asarray([self.objective_function(ind) for ind in population])
            best_idx = np.argmin(fitness)
            best_individual = population[best_idx]
            best_fitness = fitness[best_idx]
            if print_stat == True:

                print(f"=========Generation {generation}=========")
                print(f"Best Individual: {best_individual}")
                print(f"Best Score: {best_fitness}\n")

        return best_individual, best_fitness
    
    def slice_hypercube(self,lower_bounds, upper_bounds, interval):
        dim = len(lower_bounds)
        # Create a list of arrays, each containing points spaced h apart for each dimension
        points = [np.arange(lower_bounds[i], upper_bounds[i], interval[i]) for i in range(dim)]
        
        # Use meshgrid to create a grid of points in n-dimensional space
        grids = np.meshgrid(*points, indexing='ij')
        
        # Flatten and combine the grid points into a single 2D array
        grid_points = np.vstack([grid.ravel() for grid in grids]).T
        
        # Generate all vertices for smaller hypercubes within each grid cell
        offsets = np.array(list(product(*[[0, val] for val in interval])))
        res = np.array([grid_points + offset for offset in offsets])
        return res
    
    def clustering(self):
        # length of each parts in each dimension
        inc_int = (self.boundaries[:,1]-self.boundaries[:,0])/self.parts

        # Hypercubes
        hypercubes_edges = self.slice_hypercube(lower_bounds=self.boundaries[:,0],
                                                upper_bounds=self.boundaries[:,1],
                                                interval=inc_int)  
        
        cluster = []
        for hypercube_id in range(hypercubes_edges.shape[1]):
            X0 = hypercubes_edges[:,hypercube_id,:]
            F_list = self.systemeq(X0.T)

            # cek jika f yang berubah tanda dari F_list jika dievaluasi di tiap edge hypercube
            product_combination = np.array([[a*b for a,b in combinations(F_list[i],2)] for i in range (F_list.shape[0])])

            # jika semua f dari F_list berubah tanda jika dievaluasi di tiap edge hypercube, maka ada akar di situ
            change_sign = np.array([np.any(product_combination[i]<0) for i in range (product_combination.shape[0])])
            if np.all(change_sign==True):
                # print(f'Ada akar di sini: \nX0={X0}')
                cluster.append(X0)

        self.cluster = np.array(cluster)
        print(f"Number of clusters containing root: {self.cluster.shape[0]}")

    def root_elimination(self,root_archive):
        if self.dim == 1:
            list_criteria = [element for sublist in root_archive for element in sublist] #convert from 2D array into 1D array
        else:
            list_criteria = root_archive
        eligible_roots = np.array([x for x in list_criteria if (self.objective_function(x))<self.condition])
        id_duplicated_roots = []
        for i in range(len(eligible_roots)):
            for j in range (i+1,len(eligible_roots)):
                if np.linalg.norm(eligible_roots[i]-eligible_roots[j])<self.delta:
                    id_duplicated_roots.append([i,j])
        id_duplicated_roots = np.unique(id_duplicated_roots,axis=0)
        deselected_id_duplicated_roots = []
        for i in range (len(id_duplicated_roots)):
            root_a = self.objective_function(eligible_roots[id_duplicated_roots[i][0]])
            root_b = self.objective_function(eligible_roots[id_duplicated_roots[i][1]])
            if root_a<=root_b:
                id_duplicated_root = id_duplicated_roots[i][1]
            else:
                id_duplicated_root = id_duplicated_roots[i][0]
            deselected_id_duplicated_roots.append(id_duplicated_root)

        if deselected_id_duplicated_roots:
            unique_roots = np.ones(len(eligible_roots),dtype=bool)
            unique_roots[deselected_id_duplicated_roots] = False
            final_root = eligible_roots[unique_roots]
        else:
            final_root = eligible_roots
        return final_root

    def GA_evaluation(self,verbose = False, superverbose = False):
        self.archive = []
        self.score = []
        for i in range (self.cluster.shape[0]):
            subbound = np.array([[self.cluster[i,:,:][:,d].min(),self.cluster[i,:,:][:,d].max()] for d in range(self.cluster.shape[2])])
            root,root_score = self.GA(population_size=self.population_size,
                                 boundaries=subbound,
                                 max_generation=self.max_generation,
                                 mutation_rate=self.mutation_rate,
                                 seed = self.seed,
                                 print_stat=superverbose)
            self.archive.append(root)
            self.score.append(root_score)
            if verbose == True:
                print(f'\n====== Cluster {i+1} ======\n')
                print(f'Roots = {self.archive}')

        self.final_root = self.root_elimination(self.archive)
