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
% % Problem 4
% boundaries = repmat([-40, 40], dim, 1);
% % Problem 7
boundaries = [0,2;-10,10;-1,1];

hdeopt = HDE(boundaries,pop_size,parts,max_gen,mutation_factor,crossover_rate,epsilon,tau_d,seed);

[final_root,final_score,num_cluster] = hdeopt.DE_evaluation(verbose,print_stat,visual_properties)
% hdeopt.visualization2D(visual_properties)
% pause(5);
% close;


%% Exporting Statistic
clear; clc;

squared_parts = [200,400];
parts = round(sqrt(squared_parts));
max_gen = [100,250,400];

epsilon = 1e-6;
tau_d = 0.1;
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

% Define boundaries
% % Problem 1 and 5
% boundaries = repmat([-10, 10], dim, 1);
% % Problem 2
% boundaries = [-1,3;-17,4];
% % Problem 4
% boundaries = repmat([-40, 40], dim, 1);
% Problem 7
boundaries = [0,2;-10,10;-1,1];

max_iter = 110;

basename = 'hde';
us = '_';
extension = '.avi';
xlsx = '.xlsx';

disp('-start-')
failedIndices = [];
for gm = 1: size(max_gen,2)
    for pt = 1: size(parts,2)
        sheet1 = []; 
        sheet2 = [];
        sheet3 = [];
        disp("----------")
        fprintf('parts=%d\nsquared_parts=%d\nmax_gen=%d\n',parts(pt),squared_parts(pt),max_gen(gm))
        for iter=1:max_iter
            try
                fprintf('Iteration: %d\n',iter)
                seed = 'shuffle';
                rngseed = rng(seed);
                print_seed = rngseed.Seed;
                visual_properties = struct('show_visual',false, ...
                'save_visual', false, ...
                'file_name', [basename,us,num2str(squared_parts(pt)),us,num2str(max_gen(gm)),us,num2str(iter),us,num2str(rngseed.Seed),extension]);
                tic;

                hdeopt = HDE(boundaries,pop_size,parts(pt),max_gen(gm),mutation_factor,crossover_rate,epsilon,tau_d,seed);
                [final_root,final_score,num_cluster] = hdeopt.DE_evaluation(verbose,print_stat);
                num_root = size(final_root,1);
                best_score = min(final_score);

                if isempty(final_root)
                    final_root = NaN(1,dim);
                    final_score = NaN;
                end

                elapsed_time = toc;

                % 1st Sheet
                num_iter = iter*ones(size(final_root,1),1);
                sheet1 = [sheet1; num_iter,final_root,final_score'];

                % 2nd Sheet
                sheet2 = [sheet2; iter,num_root,best_score,num_cluster,elapsed_time];
                sheet3 = [sheet3; iter,print_seed];
                writematrix(sheet1 ,[basename,us,num2str(squared_parts(pt)),us,num2str(max_gen(gm)),xlsx],'Sheet','Final Root')
                writematrix(sheet2,[basename,us,num2str(squared_parts(pt)),us,num2str(max_gen(gm)),xlsx],'Sheet','Statistic')
                writematrix(sheet3,[basename,us,num2str(squared_parts(pt)),us,num2str(max_gen(gm)),xlsx],'Sheet','Random Seed')
                close;
            catch ME
                % If there is an error, display a warning and skip to the next iteration
                warning('Error untuk pt %d - gm %d - iter %d: %s', squared_parts(pt),max_gen(gm),iter, ME.message);
                % Store the index of the failed computation
                failedIndices = [failedIndices; [squared_parts(pt),max_gen(gm),iter]];
                continue;
            end
        end
    end
end
writematrix(failedIndices,['failed_indices',xlsx],'Sheet','Failed Indices')
disp('-end-')

%% Experimentation
% clc;
% he = hdeopt.slice_hypercube(boundaries(:,1),boundaries(:,2),parts)
% cluster = {};
% for hypercube_id = 1:size(he, 2)
%     X0 = he(:,hypercube_id,:);
%     X0 = reshape(X0, [], size(X0, 3));
% 
%     F_list = zeros(dim,size(X0,1));
%     for i = 1:size(X0,1)
%         F_list(:,i) = hdeopt.system_equations(X0(i,:));
%     end
% 
%     product_combination = zeros(size(F_list, 1), nchoosek(size(F_list, 2), 2));
%     for i = 1:size(F_list, 1)
%         combinations = nchoosek(F_list(i,:), 2);
%         product_combination(i,:) = prod(combinations, 2)';
%     end
%     change_sign_1 = any(product_combination < 0, 1);
%     change_sign_2 = any(product_combination < 0, 2);
%     if all(change_sign_2)
%         cluster{end+1} = X0;
%     end
% end
% cluster = cat(3, cluster{:})


% %%
% clc;
% X0 = cluster(:,:,3) %3,6,8,9,13,14
% 
% distances = squareform(pdist(X0))
% 
% dist_pair = [];
% for u=1:2^dim
%     [~,si]=sort(distances(u,:));
%     v = si(2:1+dim);
%     [U, V] = meshgrid(u, v);
%     combined_array = [U(:), V(:)];
%     dist_pair = [dist_pair;combined_array];
% end
% 
% sorted_dist_pair = sort(dist_pair, 2);
% unique_dist_pair = unique(sorted_dist_pair, 'rows')
% product_edge = zeros(dim,size(unique_dist_pair,1));
% for i = 1:size(unique_dist_pair,1)
%     fun1 = hdeopt.system_equations(X0(unique_dist_pair(i,1),:));
%     fun2 = hdeopt.system_equations(X0(unique_dist_pair(i,2),:));
%     product_edge(:,i) = fun1.*fun2;
% end
% product_edge
% product_edge<0
% a = any(product_edge<0,1)
% 
% %%
% 
% u = 1;
% v = [2, 4, 6];
% 
% [U, V] = meshgrid(u, v)
% 
% % Combine U and V into a 2D array
% combined_array = [U(:), V(:)]
% 
% %%
% clc;
% X0 = cluster(:,:,1) %3,6,8,9,13,14
% % X0 = [-10,-10;-10,-9;-9,-9;-9,-10]
% sides = getHypercubeSides(X0)
% f_value = zeros(dim,size(sides,3));
% i = 1;
% for d3 = 1:size(sides,3)
%     f_value(:,d3) = hdeopt.system_equations(sides(1,:,d3)).*hdeopt.system_equations(sides(2,:,d3));
% end
% 
% f_value
% 
% 
% %%
% clc;
% 
% ideal = [0.75,0.875;-1.25,0;-0.25,-0.125];
% 
% hpc = [
%     0.75,-1.25,-0.25;
%     0.75,-1.25,-0.125;
%     0.75,0,-0.25;
%     0.75,0,-0.125;
%     0.875,-1.25,-0.25;
%     0.875,-1.25,-0.125;
%     0.875,0,-0.25;
%     0.875,0,-0.125
%     ]
% 
% sides = getHypercubeSides(hpc);
% f_value = zeros(dim,size(sides,3));
% for d3 = 1:size(sides,3)
%     f_value(:,d3) = hdeopt.system_equations(sides(1,:,d3)).*hdeopt.system_equations(sides(2,:,d3));
% end
% f_value