%% Trial
clear; clc; close all;

m_cluster = 400;
k_max = 400;

dim = 3;

% Problem 1
dim = 2;
mutation_factor=0.864519198;
crossover_rate=0.860650007;
tau_d=0.4;
gamma = -0.2;
epsilon = 1e-6;
k_cluster = 10;
m = 250;
boundaries = repmat([-10, 10], dim, 1);
verbose = true;
print_stat = false;
seed = 'shuffle';


visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'sdde.avi');

% Define boundaries
% Problem 1 and 5
% boundaries = repmat([-10, 10], dim, 1);
% % Problem 2
% boundaries = [-1,3;-17,4];
% % Problem 7
% boundaries = [0,2;-10,10;-1,1];

CDE = CDE(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,tau_d,gamma,mutation_factor,crossover_rate,seed);

% sddeopt.initialization();
% sddeopt.cluster_center
% sddeopt.cluster_radius
% pop = sddeopt.cluster_iter_points
% sddeopt.cluster_iter_bestsol

[final_root,final_score,num_cluster] = CDE.DE_evaluation(verbose, print_stat,visual_properties)

