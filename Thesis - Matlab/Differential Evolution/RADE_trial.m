clear; clc

pop_size=100;
max_gen=100;
F_init=0.5;
CR_init=0.5;
num_l=20;
theta=1e-3;
tau_d=0.4;
s_max=100;
print_gen=true;
Hm = 100;
dim = 2;
seed = 'shuffle';
beta = 100;
rho = 1e-8;

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

radeopt = RADE(boundaries,pop_size,num_l,max_gen,s_max,theta,tau_d,F_init,CR_init,Hm,beta,rho,seed);

rng(seed);
population = radeopt.generate_points(pop_size,boundaries,seed)

fitness = zeros(1, radeopt.population_size);
for i = 1:radeopt.population_size
    fitness(i) = radeopt.objective_function(population(i, :));
end
[best_fitness, best_idx] = min(fitness);
best = population(best_idx, :);
subpop = zeros(radeopt.num_per_subpopulation,radeopt.dim,radeopt.population_size);
for i = 1:radeopt.population_size
    subpop(:, :, i) = radeopt.subpopulating(population(i, :), population, radeopt.num_per_subpopulation);
end
k=1;
for gen = 1:radeopt.max_generation
    fprintf('=====Generation %d=====\n',gen)
    S_F = [];
    S_CR = [];
    for i = 1:radeopt.population_size
        [F_i, CR_i] = radeopt.update_parameter();
        x_i = population(i, :);
        mutant = radeopt.mutation_penalty(x_i, subpop(:, :, i), radeopt.boundaries, F_i);
        trial = radeopt.crossover(population(i, :), mutant, CR_i);
        trial_fitness = radeopt.fitness_function(trial);
        fprintf("fitness(%d) = %.5f \n",i,fitness(i))
        fprintf("trial_fitness = %.5f \n",trial_fitness)
        
        if trial_fitness < fitness(i)
            fitness(i) = trial_fitness;
            population(i, :) = trial;
            archive = radeopt.update_archive(trial);
            S_F = [S_F, F_i];
            S_CR = [S_CR, CR_i];
            if trial_fitness < best_fitness
                best_idx = i;
                best = trial;
            end
        end

        if ~isempty(S_F) && ~isempty(S_CR)
            radeopt.update_history(S_F, S_CR, k);
            k = k + 1;
            if k > radeopt.max_memories_size
                k = 1;
            end
        end
    end
end

radeopt.archive