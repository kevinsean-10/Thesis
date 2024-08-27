%% Trial
clear;clc;close all;

sq_pt = 250;
parts = round(sqrt(sq_pt));
max_gen = 400;

epsilon = 1e-6;
tau_d = 0.01;
pop_size = 300;
dim = 3;
mutation_factor=0.831987414;
crossover_rate=0.887292888;
seed = 'shuffle';
print_stat = false;
verbose = false;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'hde.avi');


print_stat = false;
verbose = true;

rngseed = rng(seed);
disp(rngseed.Seed)
basename = 'hde';
us = '_';
extension = '.avi';

visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', [basename,us,num2str(sq_pt),us,num2str(max_gen),us,num2str(rngseed.Seed),extension]);


% Define boundaries
% % Problem 1 and 5
% boundaries = repmat([-10, 10], dim, 1);
% % Problem 2
% boundaries = [-1,3;-17,4];
% % Problem 3
boundaries = [0,2;-10,10;-1,1];

hdeopt = HDE(boundaries,pop_size,parts,max_gen,mutation_factor,crossover_rate,epsilon,tau_d,seed);

[final_root,final_score,num_cluster] = hdeopt.DE_evaluation(verbose,print_stat,visual_properties)
% hdeopt.visualization2D(visual_properties)
% pause(5);
% close;

