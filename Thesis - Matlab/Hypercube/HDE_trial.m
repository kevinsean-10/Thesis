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
verbose = false;
visual_properties = struct('show_visual',true, ...
    'save_visual', false, ...
    'file_name', 'hde.avi');

% how many parts/slices do you desire in each dimension?
parts = 100;

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

hdeopt = HDE(boundaries,m,parts,gen_max,mutation_factor,crossover_rate,epsilon,delta,seed);

[final_root,final_score] = hdeopt.DE_evaluation(verbose,print_stat)
hdeopt.visualization2D(visual_properties)


%% Exporting Statistic
clear; clc;

epsilon = 1e-5;
delta = 0.01;
m = 100;
gen_max = 300;
dim = 2;
mutation_factor=0.1;
crossover_rate=0.5;
print_stat = false;
verbose = false;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'hde.avi');
parts = 100;
boundaries = repmat([-10, 10], dim, 1);

max_iter = 10;

disp('-start-')
sheet1 = []; 
sheet2 = [];
for iter=1:max_iter
    fprintf('Iteration: %d\n',iter)
    tic;
    seed = 'shuffle';

    hdeopt = HDE(boundaries,m,parts,gen_max,mutation_factor,crossover_rate,epsilon,delta,seed);
    [final_root,final_score] = hdeopt.DE_evaluation(verbose,print_stat);

    elapsed_time = toc;

    % 1st Sheet
    num_iter = iter*ones(size(final_root,1),1);
    sheet1 = [sheet1; num_iter,final_root,final_score'];

    % 2nd Sheet
    num_root = size(final_root,1);
    best_score = min(final_score);
    sheet2 = [sheet2; iter,num_root,best_score,elapsed_time];
    writematrix(sheet1 ,'HDE.xlsx',Sheet='Final Root')
    writematrix(sheet2,'HDE.xlsx',Sheet='Statistic')
end

disp('-end-')
