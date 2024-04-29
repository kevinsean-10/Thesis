%% Trial
clear; clc;

dim = 2;
mutation_factor=0.1;
crossover_rate=0.4;
tau_d=0.4;
m_cluster = 250;
gamma = -0.2;
epsilon = 1e-6;
delta = 0.01;
k_cluster = 10;
m = 250;
k_max = 250;
seed = 'shuffle';
verbose = true;
print_stat = false;
visual_properties = struct('show_visual',true, ...
    'save_visual', false, ...
    'file_name', 'sdde.avi');

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

sddeopt = SDDE(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_factor,crossover_rate,seed);

% sddeopt.DE_evaluation(verbose,print_stat,visual_properties)

%% Exporting Statistic
clear; clc;

dim = 2;
mutation_factor=0.1;
crossover_rate=0.4;
tau_d=0.4;
m_cluster = 250;
gamma = -0.2;
epsilon = 1e-6;
delta = 0.01;
k_cluster = 10;
m = 250;
k_max = 250;
boundaries = repmat([-10, 10], dim, 1);
verbose = false;
print_stat = false;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'sdde.avi');

max_iter = 1;

disp('-start-')
sheet1 = []; sheet2 = [];
for iter=1:max_iter
    fprintf('Iteration: %d\n',iter)
    tic;
    seed = 'shuffle';

    sddeopt = SDDE(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_factor,crossover_rate,seed);
    [final_root,final_score] = sddeopt.DE_evaluation(verbose,print_stat,visual_properties);

    elapsed_time = toc;

    % 1st Sheet
    num_iter = iter*ones(size(final_root,1),1);
    sheet1 = [sheet1; num_iter,final_root,final_score'];

    % 2nd Sheet
    num_root = size(final_root,1);
    best_score = min(final_score);
    sheet2 = [sheet2; iter,num_root,best_score,elapsed_time];
    writematrix(sheet1 ,'SDDE.xlsx',Sheet='Final Root')
    writematrix(sheet2,'SDDE.xlsx',Sheet='Statistic')
    
end

disp('-end-')

%%
clc

tic;
for i=1:10
res = sddeopt.DE(boundaries,250,250,mutation_factor,crossover_rate,seed,true)
end
elapsed_time = toc