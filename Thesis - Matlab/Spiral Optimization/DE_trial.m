clear;clc;close all;
dim = 2;
boundaries = repmat([-10, 10], dim, 1);
max_iter = 100;
population_size = 100;
scaling_factor = 0.8;
crossover_prob = 0.8;
seed = 'shuffle';

rngseed = rng(seed);
disp(rngseed.Seed)
basename = 'hde';
us = '_';
extension = '.avi';
visual_properties = struct('show_visual',true, ...
    'save_visual', false, ...
    'file_name', [basename,us,num2str(population_size),us,num2str(max_iter),us,num2str(rngseed.Seed),extension]);

deopt = DE_class(boundaries,max_iter,population_size,scaling_factor,crossover_prob,seed);
[pop,bestsol] = deopt.generate_points(population_size,boundaries,seed);
reshape([pop.Position], dim, [])';
[pop,bestsol] = deopt.DE(pop,bestsol,boundaries,max_iter,scaling_factor,crossover_prob,true,visual_properties)
