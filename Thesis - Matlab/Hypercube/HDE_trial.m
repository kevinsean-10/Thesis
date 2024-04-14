clear;clc;

epsilon = 1e-5;
delta = 0.01;
m = 100;
gen_max = 300;
dim = 2;
mutation_factor=0.1;
crossover_rate=0.5;
seed = 'shuffle';
print_stat = false;
verbose = true;

% how many parts/slices do you desire in each dimension?
parts = 100;

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

hdeopt = HDE(boundaries,m,parts,gen_max,mutation_factor,crossover_rate,epsilon,delta,seed);

hdeopt.clustering()
archive = hdeopt.DE_evaluation(verbose,print_stat)


