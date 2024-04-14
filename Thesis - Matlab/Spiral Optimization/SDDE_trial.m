clear; clc;

dim = 2;
mutation_factor=0.1;
crossover_rate=0.5;
tau_d=0.4;
m_cluster = 100;
gamma = -0.2;
epsilon = 1e-7;
delta = 0.01;
k_cluster = 100;
m = 100;
k_max = 200;
seed = 'shuffle';
verbose = false;
print_stat = false;

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

sddeopt = SDDE(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_factor,crossover_rate,seed);

sddeopt.DE_evaluation(verbose,print_stat)

