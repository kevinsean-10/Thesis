%% Trial
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
seed = 'shuffle';
verbose = false;
visual_properties = struct('show_visual',true, ...
    'save_visual', false, ...
    'file_name', 'spiral_dynamic.avi');

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

spiropt = SD(boundaries,m_cluster,k_cluster,m,k_max,epsilon,delta,gamma,theta,r,seed);
[final_root,final_score] = spiropt.spiral_opt_evaluation(verbose,visual_properties)


%% Exporting Statistic
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
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'spiral_dynamic.avi');


verbose = false;
visualize2d = false;

boundaries = repmat([-10, 10], dim, 1);

max_iter = 10;

disp('-start-')
sheet1 = []; sheet2 = [];
for iter=1:max_iter
    fprintf('Iteration: %d\n',iter)
    tic;
    seed = 'shuffle';

    spiropt = SD(boundaries,m_cluster,k_cluster,m,k_max,epsilon,delta,gamma,theta,r,seed);
    [final_root,final_score] = spiropt.spiral_opt_evaluation(verbose,visual_properties);

    if isempty(final_root)
        final_root = NaN(1,dim);
        final_score = NaN;
    end

    elapsed_time = toc;

    % 1st Sheet
    num_iter = iter*ones(size(final_root,1),1);
    sheet1 = [sheet1; num_iter,final_root,final_score'];

    % 2nd Sheet
    num_root = size(final_root,1);
    best_score = min(final_score);
    sheet2 = [sheet2; iter,num_root,best_score,elapsed_time];
    
    writematrix(sheet1 ,'SD.xlsx',Sheet='Final Root')
    writematrix(sheet2,'SD.xlsx',Sheet='Statistic')

end

disp('-end-')
