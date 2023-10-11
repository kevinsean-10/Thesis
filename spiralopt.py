import numpy as np
import sobol_seq

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

def objective_function(x): # x:tuple n-dimension
    f = 0
    """Schwefel"""
    # for i in range (len(x)):
    #     sum_sq = 0
    #     for j in range (i+1):
    #         sum_sq += x[j]
    #     f += sum_sq**2 
    # """2^n Minima"""
    # for i in range (len(x)):
    #     f += (x[i]**4-16*x[i]**2+5*x[i])
    """Rastrigin"""
    # for i in range (len(x)):
    #     f += (x[i]**2-10*np.cos(2*np.pi*x[i])+10)
    """Problem 1 Papernya Pak Kun"""
    f1 = np.exp(x[0]-x[1])-np.sin(x[0]+x[1])
    f2 = (x[0]*x[1])**2-np.cos(x[0]+x[1])
    f = 1/(1+np.abs(f1)+np.abs(f2))
    return f

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
        for j in range (0,i+1):
            product = np.eye(dim)
            product *= generate_Rij(dim-i-1,dim+1-j-1,dim,theta)
        Rn *= product
    return Rn

def maximize(set_of_points):
    z = []
    z_max = 0
    for i in range (len(set_of_points)):
        z.append(objective_function(set_of_points.T)[i])
        if z[i]>z_max:
            z_max = z[i]
            idx_max = i
    x_max = set_of_points[idx_max]
    return z_max,idx_max,x_max

def update_point(set_of_points,Sn,dim):
    (z_star,idx_star,x_star) = maximize(set_of_points)
    new_set_of_points = np.copy(set_of_points)
    for i in range (len(new_set_of_points)):
        # perkalian matriks
        poin = np.dot(Sn,set_of_points[i].reshape(-1,1)) - np.dot((Sn-np.identity(dim)),x_star.reshape(-1,1))
        new_set_of_points[i] = poin.T
    return new_set_of_points

def iter_error(set_of_points,iter,npoint):
    err = 0
    for i in range (npoint):
        diff = np.abs(np.sqrt(np.sum(num**2 for num in set_of_points[iter][i]))-np.sqrt(np.sum(num**2 for num in set_of_points[iter-1][i])))
        if diff>err:
            err = diff
    return err
    
def SpiralOpt(low_point,high_point, dim, npoint,r = 0.95,theta=np.pi/4, iter_max=100, error_max = 10**(-5),random=0, show_err=False, show_objective_function=False):
    np.random.seed(random)
    iter_points = {}
    iter = 0
    iter_points[iter] = generate_points(dim,npoint,low_point,high_point)
    Rn = generate_Rn(dim,theta)
    Sn = r*Rn
    while iter <= iter_max :
        iter_points[iter+1] = update_point(iter_points[iter],Sn,dim)
        error = iter_error(iter_points,iter+1,npoint)
        if error < error_max:
            break
        iter += 1

    return_points = [iter_points[iter][:,i].mean() for i in range(dim)]
    if show_err == True:
        print(error)
    obj_fun = objective_function(iter_points[iter].T).mean()
    if show_objective_function ==True:
        print(obj_fun)
    return return_points, obj_fun