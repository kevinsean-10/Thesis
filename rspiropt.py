import numpy as np
import sobol_seq

# def generate_points(dim,npoint,low=-10,high=10):
#     if type(low) != type(high):
#         raise TypeError('The type of "low" and "high" should be the same.')
#     if type(low) == int:
#         boundaries = [(low,high) for _ in range (dim)]
#     elif type(low) == list or type(low) == np.ndarray:
#         if len(low) != len(high):
#             raise TypeError('The length of "low" and "high" should be the same.')
#         else:
#             boundaries = [(low[i],high[i]) for i in range (len(low))]

#     # Generate Sobol sequence points
#     sobol_points = sobol_seq.i4_sobol_generate(dim, npoint)

#     # Scale the Sobol points to fit within the specified boundaries
#     scaled_points = []
#     for i in range(dim):
#         a, b = boundaries[i]
#         scaled_dim = a + sobol_points[:, i] * (b - a)
#         scaled_points.append(scaled_dim)

#     # Trane the scaled points to get points per dimension
#     scaled_points = np.array(list(map(list, zip(*scaled_points))))
#     return scaled_points

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

# def maximize(set_of_points,objective_function):
#     z = []
#     z_max = 0
#     for i in range (len(set_of_points)):
#         z.append(objective_function(set_of_points.T)[i])
#         if z[i]>z_max:
#             z_max = z[i]
#             idx_max = i
#     x_max = set_of_points[idx_max]
#     return z_max,idx_max,x_max

def update_point(set_of_points,Sn,dim,objective_function):
    (z_star,idx_star,x_star) = maximize(set_of_points,objective_function)
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
    
def SpiralOpt(low_point,high_point,objective_function, dim, npoint,r = 0.95,theta=np.pi/4, iter_max=100, error_max = 10**(-5),random=0, show_err=False, show_objective_function=False):
    np.random.seed(random)
    iter_points = {}
    iter = 0
    iter_points[iter] = generate_points(dim,npoint,low_point,high_point)
    Rn = generate_Rn(dim,theta)
    Sn = r*Rn
    while iter <= iter_max :
        iter_points[iter+1] = update_point(iter_points[iter],Sn,dim,objective_function)
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

    # Trane the scaled points to get points per dimension
    scaled_points = np.array(list(map(list, zip(*scaled_points))))
    return scaled_points

"""MAXIMIZE FUNCTION"""
def maximize(set_of_points,objective_function):
    z = []
    z_max = 0
    F = objective_function(set_of_points.T)
    for i in range (len(set_of_points)):
        z.append(F[i])
        if z[i]>z_max:
            z_max = z[i]
            idx_max = i
    x_max = set_of_points[idx_max]
    return z_max,idx_max,x_max

"""FUNCTION CLUSTER"""
def function_cluster(y,lendict,objective_function,cluster_center,cluster_radius):
    min_dist_cluster = 10**100
    for ci,cc in cluster_center.items():
        dist = np.linalg.norm(cc-y)
        if dist<=min_dist_cluster:
            xc = cc
            cluster_id = ci
            min_dist_cluster = dist
    xt = (xc + y)/2
    # print(xt,xc,y)
    Fxt = objective_function(xt)
    Fxc = objective_function(xc)
    Fy = objective_function(y)
    # print(Fxt,Fxc,Fy)
    if Fxt < Fy and Fxt < Fxc:
        cluster_center[lendict] = y
        cluster_radius[lendict] = np.linalg.norm(y-xt)
    elif Fxt > Fy and Fxt > Fxc:
        cluster_center[lendict] = y
        cluster_radius[lendict] = np.linalg.norm(y-xt)
        function_cluster(xt,lendict+1,objective_function,cluster_center,cluster_radius)
    elif Fy > Fxc:
        cluster_center[cluster_id] = y

    cluster_radius[cluster_id] =  np.linalg.norm(y-xt)
    # return cluster_center,cluster_radius

def cluster_boundaries(center_point,radius):
    return np.array([[center_point[i]-radius,center_point[i]+radius]for i in range (len(center_point))])

def root_SpiralOpt(objective_function,m_cluster,gamma,epsilon,delta,k_cluster,m,r,theta,k_max,dim):
    k=0
    iter_points = {}
    boundaries = np.array([(-10,10) for _ in range (dim)])
    iter_points[k] = generate_points(dim,m_cluster,boundaries[:,0],boundaries[:,1])
    x_prime = maximize(iter_points[0],objective_function)[-1]
    min_boundaries = 10**100
    for i in range (len(boundaries)):
        abs_disc = np.abs(boundaries[i,1]-boundaries[i,0])
        if abs_disc<=min_boundaries:
            min_boundaries = abs_disc
    radius = min_boundaries/2
    cluster_center,cluster_radius = {},{}
    cluster_center[0],cluster_radius[0] = x_prime,radius
    Sn = r*generate_Rn(dim=dim,theta=theta)
    while k<k_cluster:
        potential_cluster_center = []
        F = objective_function(iter_points[k].T)
        for i in range (m_cluster):
            if F[i] > gamma:
                potential_cluster_center.append(iter_points[k][i])
        # len(potential_cluster_center)
        # print('')
        for i in range (len(potential_cluster_center)):
            # print(f'Titik ke-{i}')
            lendict = len(list(cluster_center.keys()))
            function_cluster(potential_cluster_center[i],lendict,objective_function,cluster_center,cluster_radius)
            # print(cluster_center,cluster_radius)
            # print('')
        iter_points[k+1]=update_point(iter_points[k],Sn,dim,objective_function=objective_function)
        k+=1
    roots = []
    roots_values = []
    for i in range (len(cluster_center)):
        bound = cluster_boundaries(cluster_center[i],cluster_radius[i])
        lp = bound[:,0]
        hp = bound[:,1]
        root,value = SpiralOpt(lp,hp,objective_function,dim,npoint=m,r = r,theta=theta, iter_max=k_max, error_max = 10**(-5),random=0, show_err=False, show_objective_function=False)
        roots.append(root)
        roots_values.append(value)
    print(roots)
    print(roots_values)
    eligible_roots = np.array([x for x in roots if (1-objective_function(x))<epsilon])
    duplicated_roots = []
    for i in range(len(eligible_roots)):
        for j in range (i+1,len(eligible_roots)):
            if np.linalg.norm(eligible_roots[i]-eligible_roots[j])<delta:
                duplicated_roots.append([eligible_roots[i],eligible_roots[j]])
    duplicated_roots = np.unique(duplicated_roots,axis=0)
    # print(duplicated_roots)
    deselected_duplicated_roots = []
    for i in range (len(duplicated_roots)):
        root_a = objective_function(duplicated_roots[i][0])
        root_b = objective_function(duplicated_roots[i][1])
        if root_a>root_b:
            deselected_duplicated_roots.append(list(duplicated_roots[i][1]))
        else:
            deselected_duplicated_roots.append(list(duplicated_roots[i][0]))
    # print(deselected_duplicated_roots)

    if deselected_duplicated_roots:
        exclude_condition = np.all(eligible_roots != np.array(deselected_duplicated_roots)[:, np.newaxis], axis=2).all(axis=0)
        final_root = eligible_roots[exclude_condition]
    else:
        final_root = eligible_roots
    return final_root
