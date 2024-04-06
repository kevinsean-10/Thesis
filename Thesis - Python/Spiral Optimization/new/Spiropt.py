import numpy as np
from scipy.stats import qmc

"""GENERATE POINTS USING SOBOL SEQUENCE"""
def generate_points(dim: int, npoint: int, low=-10, high=10, sobol=True, seed=None):
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
        boundaries = [(low, high) for _ in range(dim)]
    elif type(low) in (list, np.ndarray):
        if len(low) != len(high):
            raise TypeError('The length of "low" and "high" should be the same.')
        else:
            boundaries = [(low[i], high[i]) for i in range(len(low))]

    # Generate points based on the sobol flag
    if sobol:
        sampler = qmc.Sobol(d=dim,scramble=True,seed=seed)
        sample = sampler.random(n=npoint)
        scaled_points = qmc.scale(sample=sample,l_bounds=low,u_bounds=high)

    else:
        # Generate random points
        np.random.seed(seed)
        scaled_points = np.zeros((npoint, dim))
        for i in range(dim):
            min_val, max_val = boundaries[i]
            scaled_points[:, i] = np.random.uniform(min_val, max_val, npoint)

    return scaled_points

def generate_Rij(i,j,dim,theta):
    Rn_ij= np.eye(dim)
    Rn_ij[i-1,i-1] = np.cos(theta)
    Rn_ij[i-1,j-1] = -np.sin(theta)
    Rn_ij[j-1,i-1] = np.sin(theta)
    Rn_ij[j-1,j-1] = np.cos(theta)
    return Rn_ij

def generate_Rn(dim,theta):
    Rn = np.eye(dim)
    for i in range(0,dim):
        product = np.eye(dim)
        for j in range (0,i+1):
            product *= generate_Rij(dim-i-1,dim+1-j-1,dim,theta)
        Rn *= product
    return Rn

def update_point(set_of_points,objective_function,Sn):
    fitness = np.asarray([objective_function(ind) for ind in set_of_points])
    i_g = np.argmin(fitness)
    x_i_g = set_of_points[i_g]

    new_set_of_points = np.copy(set_of_points)
    dim = set_of_points.shape[1]

    for i in range(len(new_set_of_points)):
        poin = np.dot(Sn,set_of_points[i].reshape(-1,1)) - np.dot((Sn-np.identity(dim)),x_i_g.reshape(-1,1))
        new_set_of_points[i] = poin.T
    return new_set_of_points

def iter_error(old_set_of_points,new_set_of_points):
    err = 0
    for i in range (old_set_of_points.shape[0]):
        diff = np.abs(np.linalg.norm(old_set_of_points[i]) - np.linalg.norm(new_set_of_points[i]))
        if diff>err:
            err = diff
    return err

def spiral_opt(objective_function,boundaries,n_point,theta,radius,max_iter,max_error,seed):
    dim = boundaries.shape[0]
    Rn = generate_Rn(dim,theta)
    Sn = radius*Rn
    iter = 0
    iter_points = {}
    iter_points[iter] = generate_points(dim,n_point,boundaries[:,0],boundaries[:,1],sobol=True,seed=seed)
    while iter <= max_iter :
        iter_points[iter+1] = update_point(iter_points[iter],objective_function,Sn)
        error = iter_error(iter_points[iter],iter_points[iter+1])
        if error < max_error:
            break
        iter += 1

    return_points = np.array([iter_points[iter][:,i].mean() for i in range(dim)])
    return_points_value = objective_function(return_points)
    return return_points, return_points_value
