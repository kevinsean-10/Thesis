import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import qmc

class Spiral_Optimization:
    def __init__(
            self,
            system_equation,
            boundaries,
            m_cluster,
            k_cluster,
            m,
            k_max,
            radius,
            theta,
            gamma,
            epsilon,
            delta,
            seed = None) -> None:
        self.boundaries = boundaries
        self.systemeq = system_equation
        self.m_cluster = m_cluster
        self.k_cluster = k_cluster
        self.m = m
        self.k_max = k_max
        self.radius = radius
        self.theta = theta
        self.gamma = gamma
        self.epsilon = epsilon
        self.delta = delta
        self.seed = seed
        self.dim = boundaries.shape[0]
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
            sobol=True):
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
    
    def generate_Rij(self,i,j):
        Rn_ij= np.eye(self.dim)
        Rn_ij[i-1,i-1] = np.cos(self.theta)
        Rn_ij[i-1,j-1] = -np.sin(self.theta)
        Rn_ij[j-1,i-1] = np.sin(self.theta)
        Rn_ij[j-1,j-1] = np.cos(self.theta)
        return Rn_ij
    
    def generate_Rn(self):
        Rn = np.eye(self.dim)
        for i in range(0,self.dim):
            product = np.eye(self.dim)
            for j in range (0,i+1):
                product *= self.generate_Rij(self.dim-i-1,self.dim+1-j-1)
            Rn *= product
        return Rn
    
    def update_point(self,set_of_points,objective_function):
        Sn = self.radius*self.generate_Rn()
        fitness = np.asarray([objective_function(ind) for ind in set_of_points])
        i_g = np.argmin(fitness)
        x_i_g = set_of_points[i_g]

        new_set_of_points = np.copy(set_of_points)
        dim = set_of_points.shape[1]

        for i in range(len(new_set_of_points)):
            poin = np.dot(Sn,set_of_points[i].reshape(-1,1)) - np.dot((Sn-np.identity(dim)),x_i_g.reshape(-1,1))
            new_set_of_points[i] = poin.T
        return new_set_of_points
    
    def iter_error(self,old_set_of_points,new_set_of_points):
        err = 0
        for i in range (old_set_of_points.shape[0]):
            diff = np.abs(np.linalg.norm(old_set_of_points[i]) - np.linalg.norm(new_set_of_points[i]))
            if diff>err:
                err = diff
        return err
    
    def spiral_opt(self,
                   objective_function,
                   boundaries,
                   n_point,
                   max_iter,
                   max_error):
        dim = boundaries.shape[0]
        iter = 0
        spiral_points = {}
        spiral_points[iter] = self.generate_points(n_point,boundaries[:,0],boundaries[:,1],sobol=True)
        while iter <= max_iter :
            spiral_points[iter+1] = self.update_point(spiral_points[iter],objective_function)
            error = self.iter_error(spiral_points[iter],spiral_points[iter+1])
            if error < max_error:
                break
            iter += 1

        return_points = np.array([spiral_points[iter][:,i].mean() for i in range(dim)])
        return_points_value = objective_function(return_points)
        return return_points, return_points_value
    
    def initialization(self):
        self.iter_points = {}
        self.iter_points[0] = self.generate_points(self.m_cluster,self.boundaries[:,0],self.boundaries[:,1])

        fitness = np.asarray([self.objective_function(ind) for ind in self.iter_points[0]])
        best_idx = np.argmin(fitness)
        x_prime = self.iter_points[0][best_idx]

        radius = (self.boundaries[:,1]-self.boundaries[:,0])/2
        id_rad = np.argmin(radius)
        radii = radius[id_rad]

        self.cluster_center, self.cluster_radius = np.array([x_prime]),np.array([radii])

    """FUNCTION CLUSTER"""
    def function_cluster(self,y):
        
        dist_list = np.linalg.norm(self.cluster_center-y,axis=1)
        min_dist_id = np.argmin(dist_list)
        min_dist = dist_list[min_dist_id]
        xc = self.cluster_center[min_dist_id]
        xt = (xc + y)/2

        Fxt = self.objective_function(xt)
        Fxc = self.objective_function(xc)
        Fy = self.objective_function(y)

        if (Fxt > Fy) & (Fxt > Fxc):
            self.cluster_center = np.append(self.cluster_center,[y],axis=0)
            self.cluster_radius = np.append(self.cluster_radius, [np.linalg.norm(y-xt)],axis=0)
        elif (Fxt < Fy) & (Fxt < Fxc):
            self.cluster_center = np.append(self.cluster_center,[y],axis=0)
            self.cluster_radius = np.append(self.cluster_radius, [np.linalg.norm(y-xt)],axis=0)
            self.function_cluster(xt)
        elif Fy < Fxc:
            self.cluster_center[min_dist_id] = y
        
        self.cluster_radius[min_dist_id] =  np.linalg.norm(y-xt)

        # # update radii if the existing corresponding cluster radii is larger than tha candidate
        # if self.cluster_radius[min_dist_id] > np.linalg.norm(y-xt):
        #     self.cluster_radius[min_dist_id] =  np.linalg.norm(y-xt)
        
        # return self.cluster_center,self.cluster_radius

    def clustering(self):
        k = 0
        while k<self.k_cluster:
            potential_cluster_center = []
            F = self.objective_function(self.iter_points[k].T)
            for i in range (self.m_cluster):
                # If F(x_i)<gamma and x_i is not the center of existing cluster, x_i may have a possibility to become a cluster center
                if len(self.iter_points[k].T) == 1:
                    fungam = F[0][i]
                else:
                    fungam = F[i]
                exist_in_cluster_center = any(np.linalg.norm(self.iter_points[k][i] - ctr) < self.epsilon for ctr in self.cluster_center)
                if (fungam < self.gamma) & (exist_in_cluster_center==False):
                    potential_cluster_center.append(self.iter_points[k][i])
                
            # Apply function cluster
            for i in range (len(potential_cluster_center)):
                self.function_cluster(potential_cluster_center[i])

            self.iter_points[k+1] = self.update_point(set_of_points=self.iter_points[k],
                                                      objective_function=self.objective_function)
            k+=1
        
    def cluster_2Dvisualization(self):
        if self.dim != 2:
            print(f"Dimension {self.dim} can be visualized using cluster_visualization2D.")
        """Visualization"""
        fig, ax = plt.subplots(figsize=(10,10))
        for center,radius in zip(self.cluster_center,self.cluster_radius):
            circle = plt.Circle(center, radius, fill=False, linestyle='dotted', edgecolor='b')
            ax.add_artist(circle)

        # Set axis limits
        ax.set_xlim(self.boundaries[0])
        ax.set_ylim(self.boundaries[1])
        # ax.autoscale_view()

        # # Add labels (optional)
        # for i, center in cluster_center.items():
        #     ax.text(center[0], center[1], f'Cluster {i+1}', ha='center', va='bottom')

        # Add a title and labels (optional)
        ax.set_title('Cluster Visualization')
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')

        # Show the plot
        plt.gca().set_aspect('equal', adjustable='box')  # Make the aspect ratio equal
        plt.grid(True)
        plt.show()

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
    
    def spiral_opt_evaluation(self, verbose = False):
        self.archive = []
        self.score = []
        for i in range (len(self.cluster_center)):
            subbound = np.array([[c-self.cluster_radius[i],c+self.cluster_radius[i]] for c in self.cluster_center[i]])
            root,root_score = self.spiral_opt(objective_function=self.objective_function,
                                              boundaries=subbound,
                                              n_point=self.m,
                                              max_iter=self.k_max,
                                              max_error=self.epsilon)
            self.archive.append(root)
            self.score.append(root_score)
            if verbose == True:
                print(f'\n====== Cluster {i} ======\n')
                print(f'Roots = {self.archive}')
        self.final_root = self.root_elimination(self.archive)
        return self.final_root

