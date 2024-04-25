%% Trial
clear; clc;

dim = 2;
mutation_rate=0.1;
tau_d=0.4;
m_cluster = 250;
gamma = -0.2;
epsilon = 1e-7;
delta = 0.01;
k_cluster = 10;
m = 250;
k_max = 250;
seed = 'shuffle';
verbose = true;
print_stat = false;
visual_properties = struct('show_visual',true, ...
    'save_visual', true, ...
    'file_name', 'sdga.avi');

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

sdgaopt = SDGA(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_rate,seed);

sdgaopt.GA_evaluation(verbose,print_stat,visual_properties)

%% Exporting Statistic
clear; clc;

dim = 2;
mutation_rate=0.1;
tau_d=0.4;
m_cluster = 250;
gamma = -0.2;
epsilon = 1e-7;
delta = 0.01;
k_cluster = 10;
m = 250;
k_max = 250;
seed = 'shuffle';
verbose = false;
print_stat = false;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'sdga.avi');
boundaries = repmat([-10, 10], dim, 1);

max_iter = 10;

disp('-start-')
sheet1 = []; 
sheet2 = [];
for iter=1:max_iter
    fprintf('Iteration: %d\n',iter)
    tic;
    seed = 'shuffle';

    sdgaopt = SDGA(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_rate,seed);
    [final_root,final_score] = sdgaopt.GA_evaluation(verbose,print_stat,visual_properties);

    elapsed_time = toc;

    % 1st Sheet
    num_iter = iter*ones(size(final_root,1),1);
    sheet1 = [sheet1; num_iter,final_root,final_score'];

    % 2nd Sheet
    num_root = size(final_root,1);
    best_score = min(final_score);
    sheet2 = [sheet2; iter,num_root,best_score,elapsed_time];
    writematrix(sheet1 ,'SDGA.xlsx',Sheet='Final Root')
    writematrix(sheet2,'SDGA.xlsx',Sheet='Statistic')
end

disp('-end-')