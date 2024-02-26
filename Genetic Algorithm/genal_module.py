from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.termination.default import DefaultSingleObjectiveTermination
from pymoo.operators.mutation.pm import PolynomialMutation
import numpy as np
from itertools import product,combinations
import matplotlib.pyplot as plt
from typing import Callable,Tuple
from numpy.typing import NDArray
from functools import partial

ObjectiveFunctionType = Callable[[np.ndarray], np.ndarray]
BoundariesType = NDArray[np.int_]

# def root_objective_function(x,objective_function:ObjectiveFunctionType):
#     F_array = objective_function(x)
#     denom = 0
#     for f in F_array:
#         denom +=np.abs(f)
#     F = 1/(1+denom)
#     return F

def slice_hypercube(lower_bounds, upper_bounds, interval, dim):
    # Create a list of arrays, each containing points spaced h apart for each dimension
    points = [np.arange(lower_bounds[i], upper_bounds[i], interval) for i in range(dim)]
    
    # Use meshgrid to create a grid of points in n-dimensional space
    grids = np.meshgrid(*points, indexing='ij')
    
    # Flatten and combine the grid points into a single 2D array
    grid_points = np.vstack([grid.ravel() for grid in grids]).T
    
    # Generate all vertices for smaller hypercubes within each grid cell
    return np.array([grid_points + offset for offset in product([0, interval], repeat=dim)])

def GenAl(n_var,n_obj,xl,xu,objective_function,algorithmGA,termination,seed=1,verbose=False):
    class myProblem(ElementwiseProblem):
        def __init__(self):
            super().__init__(n_var=n_var,n_obj=n_obj,xl=xl,xu=xu)
        def _evaluate(self, x, out, *args, **kwargs):
            out['F'] = objective_function(x)

    problem = myProblem()
    result = minimize(problem=problem,
                    algorithm=algorithmGA,
                    termination=termination,
                    seed=seed,
                    verbose=verbose)
    return [result.X,result.F]

def root_GenAl(objective_function:ObjectiveFunctionType,
               root_objective_function,
               boundaries: BoundariesType,
               dim:int,
               gen_max:int,
               n_points:int,
               epsilon = 10**(-3),
               delta = 0.01,
               p_mutation:float=0.1,
               parts:int=100,
               seed:int=0,
               print_cluster=False
               ):
    # length of each parts in each dimension
    inc_int = (boundaries[0,1]-boundaries[0,0])/parts

    # Hypercubes
    hypercubes_edges = slice_hypercube(lower_bounds=boundaries[:,0],
                        upper_bounds=boundaries[:,1],
                        interval=inc_int,
                        dim=dim)
    
    """Clustering Technique"""
    cluster = []
    for hypercube_id in range(hypercubes_edges.shape[1]):
        X0 = hypercubes_edges[:,hypercube_id,:]
        F_list = objective_function(X0.T)

        # cek jika f yang berubah tanda dari F_list jika dievaluasi di tiap edge hypercube
        product_combination = np.array([[a*b for a,b in combinations(F_list[i],2)] for i in range (F_list.shape[0])])
        # jika semua f dari F_list berubah tanda jika dievaluasi di tiap edge hypercube, maka ada akar di situ
        change_sign = np.array([np.any(product_combination[i]<0) for i in range (product_combination.shape[0])])
        if np.all(change_sign==True):
            # print(f'Ada akar di sini: \nX0={X0}')
            cluster.append(X0)
    cluster = np.array(cluster)

    """Genetic Algorithm"""
    # GA parameter
    mutation = PolynomialMutation(prob=p_mutation)
    algorithmGAp1 = GA(pop_size=n_points,eliminate_duplicates=True,mutation=mutation)

    # termination variable
    terminationp1 = DefaultSingleObjectiveTermination(
        xtol=1e-8,
        cvtol=1e-6,
        ftol=epsilon,
        period=20,
        n_max_gen=gen_max,
        n_max_evals=100000
    )
    roots = []
    values = []
    for i in range (cluster.shape[0]):
        root,value = GenAl(2,
                           1,
                           xl=cluster[i,:,:][0],
                           xu=cluster[i,:,:][-1],
                           objective_function=root_objective_function,
                           algorithmGA=algorithmGAp1,
                           termination=terminationp1,
                           seed=seed)
        roots.append(root)
        values.append(value)
    roots = np.array(roots)
    values = np.array(values)
    if print_cluster==True:
        print(f'Number of Clusters containing root: {cluster.shape[0]}\n')
        print(f'Roots:\n{roots}\n\nValues: \n{values}')

    """Choosing Best Solution"""
    if dim == 1:
        list_criteria = [element for sublist in roots for element in sublist] #convert from 2D array into 1D array
    else:
        list_criteria = roots
    eligible_roots = np.array([x for x in list_criteria if root_objective_function(x)<epsilon])
    duplicated_roots = []
    for i in range(len(eligible_roots)):
        for j in range (i+1,len(eligible_roots)):
            if np.linalg.norm(eligible_roots[i]-eligible_roots[j])<delta:
                duplicated_roots.append([eligible_roots[i],eligible_roots[j]])
    duplicated_roots = np.unique(duplicated_roots,axis=0)

    deselected_duplicated_roots = []
    for i in range (len(duplicated_roots)):
        value_root_a = root_objective_function(duplicated_roots[i][0])
        value_root_b = root_objective_function(duplicated_roots[i][1])
        if dim == 1:
            if value_root_a<value_root_b:
                duplicated_root = duplicated_roots[i][1]
            else:
                duplicated_root = duplicated_roots[i][0]
        else:
            if value_root_a<value_root_b:
                duplicated_root = list(duplicated_roots[i][1])
            else:
                duplicated_root = list(duplicated_roots[i][0])
        deselected_duplicated_roots.append(duplicated_root)

    if dim == 1:
        # Reshape the 1D array to have one column
        deselected_duplicated_roots = np.array(deselected_duplicated_roots).reshape(-1, 1)

        # Compare the 2D array with the reshaped 1D array
        exclude_condition = np.all(eligible_roots != deselected_duplicated_roots, axis=0)

        # Use the boolean mask to filter eligible_roots
        final_root = eligible_roots[exclude_condition]
    else:
        if deselected_duplicated_roots:
            exclude_condition = np.all(eligible_roots != np.array(deselected_duplicated_roots)[:, np.newaxis], axis=2).all(axis=0)
            final_root = eligible_roots[exclude_condition]
        else:
            final_root = eligible_roots
    return final_root
