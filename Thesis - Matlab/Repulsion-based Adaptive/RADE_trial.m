%% Trial
clear; clc

pop_size=250;
max_gen=250;
F_init=0.5;
CR_init=0.5;
num_l=10;
theta=1e-7;
tau_d=0.4;
s_max=20;
print_gen=true;
Hm = 50;
dim = 2;
seed = 'shuffle';
beta = 1;
rho = 0.3;
visual_properties = struct('show_visual',true, ...
                            'save_visual', false, ...
                            'file_name', 'rade3d.avi');

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

radeopt = RADE(boundaries,pop_size,num_l,max_gen,s_max,theta,tau_d,F_init,CR_init,Hm,beta,rho,seed);

[final_root,final_score] = radeopt.DE_evaluation(print_gen,visual_properties)




%% Exporting Statistic
clear; clc;

pop_size=250;
max_gen=250;
F_init=0.5;
CR_init=0.5;
num_l=20;
theta=1e-3;
tau_d=0.4;
s_max=20;
print_gen=false;
Hm = 50;
dim = 2;
beta = 10;
rho = 1;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'rade.avi');
boundaries = repmat([-10, 10], dim, 1);

max_iter = 10;

disp('-start-')
sheet1 = []; 
sheet2 = [];
for iter=1:max_iter
    fprintf('Iteration: %d\n',iter)
    tic;
    seed = 'shuffle';

    radeopt = RADE(boundaries,pop_size,num_l,max_gen,s_max,theta,tau_d,F_init,CR_init,Hm,beta,rho,seed);
    [final_root,final_score] = radeopt.DE_evaluation(print_gen,visual_properties);

    elapsed_time = toc;

    % 1st Sheet
    num_iter = iter*ones(size(final_root,1),1);
    sheet1 = [sheet1; num_iter,final_root,final_score'];

    % 2nd Sheet
    num_root = size(final_root,1);
    best_score = min(final_score);
    sheet2 = [sheet2; iter,num_root,best_score,elapsed_time];
    writematrix(sheet1 ,'RADE.xlsx',Sheet='Final Root')
    writematrix(sheet2,'RADE.xlsx',Sheet='Statistic')
end

disp('-end-')

%%
% rng(seed);
% population = radeopt.generate_points(pop_size,boundaries,seed)
% 
% fitness = zeros(1, radeopt.population_size);
% for i = 1:radeopt.population_size
%     fitness(i) = radeopt.objective_function(population(i, :));
% end
% [best_fitness, best_idx] = min(fitness);
% best = population(best_idx, :);
% subpop = zeros(radeopt.num_per_subpopulation,radeopt.dim,radeopt.population_size);
% for i = 1:radeopt.population_size
%     subpop(:, :, i) = radeopt.subpopulating(population(i, :), population, radeopt.num_per_subpopulation);
% end
% k=1;
% for gen = 1:radeopt.max_generation
%     fprintf('=====Generation %d=====\n',gen)
%     S_F = [];
%     S_CR = [];
%     for i = 1:radeopt.population_size
%         [F_i, CR_i] = radeopt.update_parameter();
%         x_i = population(i, :);
%         mutant = radeopt.mutation_penalty(x_i, subpop(:, :, i), radeopt.boundaries, F_i);
%         trial = radeopt.crossover(population(i, :), mutant, CR_i);
%         trial_fitness = radeopt.fitness_function(trial);
%         fprintf("fitness(%d) = %.5f \n",i,fitness(i))
%         fprintf("trial_fitness = %.5f \n",trial_fitness)
% 
%         if trial_fitness < fitness(i)
%             fitness(i) = trial_fitness;
%             population(i, :) = trial;
%             archive = radeopt.update_archive(trial);
%             S_F = [S_F, F_i];
%             S_CR = [S_CR, CR_i];
%             if trial_fitness < best_fitness
%                 best_idx = i;
%                 best = trial;
%             end
%         end
% 
%         if ~isempty(S_F) && ~isempty(S_CR)
%             radeopt.update_history(S_F, S_CR, k);
%             k = k + 1;
%             if k > radeopt.max_memories_size
%                 k = 1;
%             end
%         end
%     end
% end

% radeopt.archive

%%

[X_surf, Y_surf] = meshgrid(linspace(boundaries(1,1),boundaries(1,2),50),linspace(boundaries(2,1),boundaries(2,2),50));

XY_surf = [X_surf(:),Y_surf(:)];

Z_surf = zeros(size(XY_surf,1),1);
for i=1:size(XY_surf,1)
    Z_surf(i) = radeopt.repulsion_function(XY_surf(i,:));
end

Z_surf = reshape(Z_surf, size(X_surf))

figure;
surf(X_surf, Y_surf, Z_surf);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Surface Plot of f(x, y)');
%%
% Define the range of x and y values
x_range = linspace(boundaries(1, 1), boundaries(1, 2), 50);
y_range = linspace(boundaries(2, 1), boundaries(2, 2), 50);

% Create a meshgrid of x and y values
[X_surf, Y_surf] = meshgrid(x_range, y_range);

% Flatten the meshgrid into a single vector of XY values
XY_surf = [X_surf(:), Y_surf(:)];

% Preallocate memory for Z_surf
Z_surf = zeros(size(XY_surf, 1), 1);

% Evaluate the function for each XY pair
for i = 1:size(XY_surf, 1)
    Z_surf(i) = radeopt.objective_function(XY_surf(i, :));
end

% Reshape Z_surf to match the size of X_surf and Y_surf
Z_surf = reshape(Z_surf, size(X_surf));

% Create the figure for the animation
figure;

% Create the first subplot for the surface plot
subplot(2, 1, 1);
surf(X_surf, Y_surf, Z_surf);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Surface Plot of f(x, y)');

% Set axis limits
axis tight;

% Create the second subplot for another type of plot (e.g., a scatter plot)
subplot(2, 1, 2);
scatter(X_surf(:), Y_surf(:), 10, Z_surf(:), 'filled');
xlabel('X');
ylabel('Y');
title('Scatter Plot of f(x, y)');

% Set axis limits
axis tight;

% Set the animation properties
num_frames = 100; % Number of frames in the animation
pause_time = 0.1; % Pause time between frames in seconds

% Loop for animating the surface plot
for frame = 1:num_frames
    % Update the Z data for the surface plot
    Z_surf = Z_surf + randn(size(Z_surf))*0.1; % Example: adding random noise to the surface
    
    % Update the surface plot with the new Z data
    subplot(2, 1, 1);
    surf(X_surf, Y_surf, Z_surf);
    
    % Update the scatter plot with the new Z data
    subplot(2, 1, 2);
    scatter(X_surf(:), Y_surf(:), 10, Z_surf(:), 'filled');
    
    % Pause for a short duration to create animation effect
    pause(pause_time);
end
