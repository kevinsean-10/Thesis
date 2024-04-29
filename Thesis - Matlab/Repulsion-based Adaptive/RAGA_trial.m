clear; clc

pop_size=1000;
max_gen=250;
F_init=0.5;
CR_init=0.5;
num_l=10;
theta=1e-3;
tau_d=0.4;
s_max=50;
print_gen=true;
Hm = 100;
dim = 2;
seed = 'shuffle';
beta = 1;
rho = 0.3;
visual_properties = struct('show_visual',true, ...
    'save_visual', false, ...
    'file_name', 'raga.avi');

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

ragaopt = RAGA(boundaries,pop_size,num_l,max_gen,s_max,theta,tau_d,F_init,CR_init,Hm,beta,rho,seed);
[final_root,final_score] = ragaopt.GA_evaluation(print_gen,visual_properties)

%% Exporting Statistic
clear; clc;

pop_size=250;
max_gen=250;
F_init=0.5;
CR_init=0.5;
num_l=20;
theta=1e-3;
tau_d=0.4;
s_max=50;
print_gen=false;
Hm = 50;
dim = 2;
beta = 100;
rho = 1e-8;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'raga.avi');

boundaries = repmat([-10, 10], dim, 1);

max_iter = 10;

disp('-start-')
sheet1 = []; 
sheet2 = [];
for iter=1:max_iter
    fprintf('Iteration: %d\n',iter)
    tic;
    seed = 'shuffle';

    ragaopt = RAGA(boundaries,pop_size,num_l,max_gen,s_max,theta,tau_d,F_init,CR_init,Hm,beta,rho,seed);
    [final_root,final_score] = ragaopt.GA_evaluation(print_gen,visual_properties)

    elapsed_time = toc;

    % 1st Sheet
    num_iter = iter*ones(size(final_root,1),1);
    sheet1 = [sheet1; num_iter,final_root,final_score'];

    % 2nd Sheet
    num_root = size(final_root,1);
    best_score = min(final_score);
    sheet2 = [sheet2; iter,num_root,best_score,elapsed_time];
    writematrix(sheet1 ,'RAGA.xlsx',Sheet='Final Root')
    writematrix(sheet2,'RAGA.xlsx',Sheet='Statistic')
end

disp('-end-')


%%
% rng(ragaopt.seed);
% population = ragaopt.generate_points(ragaopt.population_size,ragaopt.boundaries,ragaopt.seed);
% 
% fitness = zeros(1, ragaopt.population_size);
% for i = 1:ragaopt.population_size
%     fitness(i) = ragaopt.objective_function(population(i, :));
% end
% [best_fitness, best_idx] = min(fitness);
% best = population(best_idx, :);
% subpopulation = zeros(ragaopt.num_per_subpopulation,ragaopt.dim,ragaopt.population_size);
% for i = 1:ragaopt.population_size
%     subpopulation(:, :, i) = ragaopt.subpopulating(population(i, :), population, ragaopt.num_per_subpopulation);
% end
% 
% memory_id=1;
% 
% for gen = 1:ragaopt.max_generation 
%     S_F = [];
%     S_CR = [];
%     for i = 1:2:ragaopt.population_size
%         [F_i, CR_i] = ragaopt.update_parameter();
%         x_i = population(i, :);
%         parent = zeros(2, ragaopt.dim); % because of one-point crossover
% 
%         for j = 0:size(parent, 1)-1
%             fitness_subpopulation = zeros(1,size(subpopulation,1));
%             for k = 1:size(subpopulation,1)
%                 fitness_subpopulation(k) = ragaopt.objective_function(subpopulation(k,:,i+j));
%             end
%             [selected_subpopulation,selected_indices] = ragaopt.selection(subpopulation(:,:,i+j),fitness_subpopulation);
%             [~,min_idx_parent] = min(fitness_subpopulation(selected_indices));
%             parent(j+1,:) = selected_subpopulation(min_idx_parent,:);
%         end
%         offspring_set = ragaopt.crossover(parent(1,:),parent(2,:));
% 
%         for k=1:size(offspring_set,1)
%             trial = ragaopt.mutate(offspring_set(k,:),F_i,ragaopt.boundaries);
%             trial_fitness = ragaopt.fitness_function(trial);
%             [id_closest_trial,closest_trial] = ragaopt.closest_solution(trial,population);
%             closest_trial_fitness = ragaopt.fitness_function(closest_trial);
%             if trial_fitness < closest_trial_fitness
%                 fitness(id_closest_trial) = trial_fitness;
%                 population(id_closest_trial,:) = trial;
%                 ragaopt.archive = ragaopt.update_archive(trial);
%                 S_F = [S_F, F_i];
%                 S_CR = [S_CR, CR_i];
%                 if trial_fitness < best_fitness
%                     best_idx = i;
%                     best = trial
%                 end
%             end
%         end
%     end
%     if print_gen
%         fprintf("=========Generation %d=========\n", gen);
%         disp("Archive:");
%         disp(ragaopt.archive);
%     end
%     if ~isempty(S_F) && ~isempty(S_CR)
%         ragaopt.update_history(S_F, S_CR, memory_id);
%         memory_id = memory_id + 1;
%         if memory_id > ragaopt.max_memories_size
%             memory_id = 1;
%         end
%     end
% end

%%
clc
population_size = 10;
num_per_subpopulation = 4;
population = ragaopt.generate_points(population_size,boundaries,'shuffle')
for i = 1:2:population_size
    fprintf('i=%d\n',i)
    [F_i, CR_i] = ragaopt.update_parameter();
    x_i = population(i, :);
    parent = zeros(2, dim); % because of one-point crossover

    for j = 0:size(parent, 1)-1
        fprintf('j=%d\n',j)
        subpopulation_i = ragaopt.subpopulating(population(i+j, :), population, num_per_subpopulation);
        fitness_subpopulation_i = zeros(1,num_per_subpopulation);
        for k = 1:size(subpopulation_i,1)
            fitness_subpopulation_i(k) = ragaopt.objective_function(subpopulation_i(k,:));
        end
        [selected_subpopulation,selected_indices] = ragaopt.selection(subpopulation_i,fitness_subpopulation_i);
        [~,min_idx_parent] = min(fitness_subpopulation_i(selected_indices));
        parent(j+1,:) = selected_subpopulation(min_idx_parent,:);
    end
    parent
end

%%
clc
% Example vector
v = -rand(1,3)

% Calculate the mean and standard deviation of the vector
mean_v = mean(v);
std_v = std(v);

% Normalize the vector
normalized_v = (v - mean_v) / std_v
probability_distribution = normcdf(-normalized_v, mean_v, std_v)


%%
clc
% Example normalized vector
normalized_v = [-0.5, 0.2, -0.8];

% Define parameters for the normal distribution
mu = 0; % Mean of the distribution
sigma = 1; % Standard deviation of the distribution

% Compute the cumulative distribution function (CDF) of the normal distribution
% Note: Use the negative of the normalized value because we want higher probabilities for more negative values
probability_distribution = normcdf(-normalized_v, mu, sigma);

% Display the resulting probability distribution
disp('Probability distribution:');
disp(probability_distribution);

