clear; clc;

dim = 2;
mutation_rate=0.1;
tau_d=0.4;
m_cluster = 100;
gamma = -0.2;
epsilon = 1e-7;
delta = 0.01;
k_cluster = 200;
m = 100;
k_max = 200;
seed = 'shuffle';
verbose = false;
print_stat = false;

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

sdgaopt = SDGA(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_rate,seed);

sdgaopt.GA_evaluation(verbose,print_stat)


