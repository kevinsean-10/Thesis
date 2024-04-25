%% Trial

clear;clc

epsilon = 1e-3;
delta = 0.01;
m = 100;
gen_max = 1000;
dim = 2;
p_mutation = 0.1;
seed = 'shuffle';
print_stat = false;
verbose = false;
visual_properties = struct('show_visual',true, ...
    'save_visual', false, ...
    'file_name', 'hga.avi');

% how many parts/slices do you desire in each dimension? (even number only)
parts = 100;

boundaries = repmat([-10, 10], dim, 1);

hgaopt = HGA(boundaries,m,parts,gen_max,p_mutation,epsilon,delta,seed);

[final_root,final_score] = hgaopt.GA_evaluation(verbose,print_stat)
hgaopt.visualization2D(visual_properties)


%% Exporting Statistic
clear; clc;

epsilon = 1e-3;
delta = 0.01;
m = 100;
gen_max = 1000;
dim = 2;
p_mutation = 0.1;
seed = 'shuffle';
print_stat = false;
verbose = false;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'hga.avi');
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

    hgaopt = HGA(boundaries,m,parts,gen_max,p_mutation,epsilon,delta,seed);
    [final_root,final_score] = hgaopt.GA_evaluation(verbose,print_stat);

    elapsed_time = toc;

    % 1st Sheet
    num_iter = iter*ones(size(final_root,1),1);
    sheet1 = [sheet1; num_iter,final_root,final_score'];

    % 2nd Sheet
    num_root = size(final_root,1);
    best_score = min(final_score);
    sheet2 = [sheet2; iter,num_root,best_score,elapsed_time];
    writematrix(sheet1 ,'HGA.xlsx',Sheet='Final Root')
    writematrix(sheet2,'HGA.xlsx',Sheet='Statistic')
end

disp('-end-')





%%
% 
% clear;clc;
% % lower_bounds = boundaries(:,1); upper_bounds = boundaries(:,2); 
% % interval = (upper_bounds-lower_bounds)/(parts);
% % 
% % hypercubes_edges = hgaopt.slice_hypercube(lower_bounds,upper_bounds,parts);
% % 
% % cluster = {};
% % 
% % for hypercube_id = 1:size(hypercubes_edges, 2)
% %     X0 = hypercubes_edges(:,hypercube_id,:);
% %     X0 = reshape(X0, [], size(X0, 3));
% % 
% %     F_list = zeros(dim,size(X0,1));
% %     for i = 1:size(X0,1)
% %         F_list(:,i) = hgaopt.system_equations(X0(i,:));
% %     end
% % 
% %     product_combination = zeros(size(F_list, 1), nchoosek(size(F_list, 2), 2));
% %     for i = 1:size(F_list, 1)
% %         combinations = nchoosek(F_list(i,:), 2);
% %         product_combination(i,:) = prod(combinations, 2)';
% %     end
% % 
% %     change_sign = any(product_combination < 0, 2);
% %     if all(change_sign)
% %         cluster{end+1} = X0;
% %     end
% % end
% % % Convert cluster to array
% % cluster = cat(3, cluster{:})
% % fprintf("Number of clusters containing root: %d\n", size(cluster, 3));
% 
% % Given cell array of arrays
% cell_array = {[-6.1385, -0.1594], [-6.4349, 0.1538], [0.1651, 6.1198]};
% 
% % Convert cell array to matrix
% matrix = zeros(numel(cell_array), numel(cell_array{1}));
% for i = 1:numel(cell_array)
%     matrix(i, :) = cell_array{i};
% end
% 
% disp(matrix)
% 
% 
