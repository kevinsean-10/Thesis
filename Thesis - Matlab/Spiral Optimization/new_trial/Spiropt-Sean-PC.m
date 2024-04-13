clear
clc


m_cluster = 250;
gamma = -0.2;
epsilon = 1e-7;
delta = 0.01;
k_cluster = 10;
m = 250;
r = 0.95;
theta = pi/4;
k_max = 250;
dim = 2;
seed = 0; %'shuffle';
verbose = true;

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

spiropt = spiral_optimization(boundaries,m_cluster,k_cluster,m,k_max,epsilon,delta,gamma,theta,r,seed);

spiropt.spiral_opt_evaluation(verbose)


