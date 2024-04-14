clear; clc

pop_size=1000;
max_gen=100;
F_init=0.5;
CR_init=0.5;
num_l=5;
theta=1e-3;
tau_d=0.4;
s_max=100;
print_gen=true;
Hm = 50;
dim = 2;
seed = 'shuffle';
beta = 100;
rho = 1e-8;

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

ragaopt = RAGA(boundaries,pop_size,num_l,max_gen,s_max,theta,tau_d,F_init,CR_init,Hm,beta,rho,seed);
ragaopt.DE_evaluation(print_gen)
ragaopt.archive

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
