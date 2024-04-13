clear
clc


m_cluster = 100;
gamma = -0.2;
epsilon = 10^(-3);
delta = 0.01;
k_cluster = 100;
m = 500;
r = 0.95;
theta = pi/4;
k_max = 250;
dim = 2;

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);


sporoot = Spiral_Optimization(boundaries, m_cluster, k_cluster, m, k_max, r, theta, gamma, epsilon, delta);
% sporoot.generate_points(m_cluster,boundaries(:,1),boundaries(:,2))

sporoot.spiral_opt_evaluation()