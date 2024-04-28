classdef RAGA < handle
    properties
        boundaries
        dim
        population_size
        num_per_subpopulation
        max_generation
        theta
        tau_d
        beta
        rho
        archive
        criterion
        max_archive_size
        memories_F
        memories_CR
        max_memories_size
        seed
    end
    
    methods
        function obj = RAGA(boundaries,population_size, ...
                num_per_subpopulation,max_generation, ...
                max_archive_size,theta,tau_d,F_init,CR_init, ...
                max_memories_size,beta,rho,seed)
            obj.boundaries = boundaries;
            obj.dim = size(obj.boundaries,1);
            obj.population_size = population_size;
            obj.num_per_subpopulation = num_per_subpopulation;
            obj.max_archive_size = max_archive_size;
            obj.max_generation = max_generation;
            obj.theta = theta;
            obj.criterion = theta;
            obj.beta = beta;
            obj.rho = rho;
            obj.tau_d = tau_d;
            obj.max_memories_size = max_memories_size;
            obj.memories_F = ones(1, max_memories_size) * F_init;
            obj.memories_CR = ones(1, max_memories_size) * CR_init;
            obj.archive = [];
            obj.seed = seed;
        end
        
        function F_array = system_equations(obj,x)
            f1 = exp(x(1)-x(2)) - sin(x(1)+x(2));
            f2 = (x(1)*x(2))^2 - cos(x(1)+x(2));
            F_array = [f1; f2];
        end

        function res = objective_function(obj, x)
            F_array = obj.system_equations(x);
            res = sum(abs(F_array));
            res = -1 / (1 + res);
            obj.criterion = -1 + obj.theta;
        end

        % function res = objective_function(obj,x)
        %     res = 0;
        %     F_array = obj.system_equations(x);
        %     for i = 1:numel(F_array)
        %         res = res + (F_array(i))^2;
        %     end
        % end

        function result = chi_p(obj, delta_j, rho)
            % Characteristic Function
            if delta_j <= rho
                result = 1;
            else
                result = 0;
            end
        end

        function Rx = repulsion_function(obj, x)
            % Repulsion Function
            f_x = obj.objective_function(x);
            Rx = 0;
            for i = 1:size(obj.archive,1)
                x_star = obj.archive(i,:);
                delta_j = norm(x - x_star);
                Rx = Rx + exp(-delta_j) * obj.chi_p(delta_j, obj.rho);
            end
            Rx = Rx * obj.beta + f_x;
        end

        function f_x = fitness_function(obj, x)
            % Fitness Function
            if isempty(obj.archive)
                f_x = obj.objective_function(x);
            else
                f_x = obj.repulsion_function(x);
            end
        end

        function points = generate_points(obj,npoint,boundaries,seed)
            rng(seed)
            dimension = size(boundaries,1);
            p = sobolset(dimension);
            p = scramble(p,'MatousekAffineOwen');
            A = net(p,npoint);
            points = zeros(npoint,dimension);
            for i=1:dimension
               points(:,i)=round((boundaries(i,1)+(boundaries(i,2)-boundaries(i,1)).*A(:,i))*100)/100;
            end
        end

        function [selected_population,selected_indices] = selection(obj, population, fitness)
            pop_size = size(population, 1);
            selection_probs = 1 ./ (fitness+1); % Make sure no division by zero
            total_probs = sum(selection_probs);
            selection_probs = selection_probs / total_probs;
            selected_indices = randsample(1:pop_size, pop_size, true, selection_probs);
            selected_population = population(selected_indices, :);
        end

        function offspring_set = crossover(obj, parent1, parent2)
            dimension = length(parent1);
            crossover_point = randi([1, dimension], 1);
            offspring1 = [parent1(1:crossover_point), parent2(crossover_point+1:end)];
            offspring2 = [parent2(1:crossover_point), parent1(crossover_point+1:end)];
            offspring_set = [offspring1; offspring2];
        end

        function individual = mutate(obj, individual, mutation_rate, boundaries)
            for j = 1:length(individual)
                if rand() < mutation_rate
                    individual(j) = rand() * (boundaries(j, 2) - boundaries(j, 1)) + boundaries(j, 1);
                end
            end
        end

        function offspring_population = recombination(obj, population, mutation_rate, boundaries)
            offspring_population = [];
            pop_size = size(population, 1);
            for i = 1:2:pop_size
                parent1 = population(i, :);
                parent2 = population(i + 1, :);
                [offspring1, offspring2] = obj.crossover(parent1, parent2);
                offspring1 = obj.mutate(offspring1, mutation_rate, boundaries);
                offspring2 = obj.mutate(offspring2, mutation_rate, boundaries);
                offspring_population = [offspring_population; offspring1; offspring2];
            end
        end

        function [id_min_dist, closest_point] = closest_solution(obj,initial_point, set_of_points)
            diff = set_of_points - initial_point;
            distances = vecnorm(diff, 2, 2);
            [~, id_min_dist] = min(distances);
            closest_point = set_of_points(id_min_dist,:);
        end

        function subpopulation = subpopulating(obj, individual, population, t)
            % Calculate Euclidean distances and select t closest individuals.
            % Inputs:
            % individual: individual target
            % population: population where the target is
            % t: max number of units in a subpopulation
            % return_index: optional, flag to return indices of closest individuals
            % show_distances: optional, flag to print distances
            % Outputs:
            % closest_indices: array of closest index of the target individual
            % subpop: array of closest individuals of the target individual
        
            % Calculate the Euclidean distances from the individual to all others in the population
            distances = sqrt(sum((population - individual) .^ 2, 2));
        
            % Get the indices of the individuals with the smallest distances
            [~, sorted_indices] = sort(distances);
            closest_indices = sorted_indices(1:t);
        
            % Form the subpopulation with the closest individuals
            subpopulation = population(closest_indices, :);
        end

        function archive = update_archive(obj, x, archive)
            % Update the archive with a new individual x
            % Inputs:
            % x: individual
            % Outputs:
            % archive: updated archive
        
            f_x = obj.objective_function(x);
            s = size(archive,1); % current archive size

            if f_x < obj.criterion % x is a root
                if s == 0 % archive is empty
                    archive = [archive; x];
                    s = s + 1;
                else
                    % Find the closest solution x_prime in the archive to x in the decision space
                    dist_list = vecnorm(x - archive, 2, 2);
                    [dist_min, idx_min] = min(dist_list);
                    x_prime = archive(idx_min, :);
                    f_x_prime = obj.objective_function(x_prime);
        
                    if dist_min < obj.tau_d % x and x_prime are too close
                        if f_x < f_x_prime
                            x_prime = x;
                            archive(idx_min, :) = x_prime;
                        end
                    else
                        if s < obj.max_archive_size
                            archive = [archive; x];
                            s = s + 1;
                        else % archive is full
                            if f_x < f_x_prime
                                x_prime = x;
                                archive(idx_min, :) = x_prime;
                            end
                        end
                    end
                end
            end
        end

        function [Fi, CRi] = update_parameter(obj)
            % OUTPUT
            % Fi: scaling factor
            % CRi: crossover rate

            % Randomly select an index
            hi = randi(obj.max_memories_size);
            % Generate Fi using the Cauchy distribution with the location parameter MF(hi) and scale 0.1
            % Fi = tan(pi * (rand - 0.5)) + obj.memories_F(hi)
            Fi = obj.memories_F(hi) + 0.1*trnd(1,1);
            % Generate CRi using the Gaussian distribution with mean MCR(hi) and standard deviation 0.1
            CRi = normrnd(obj.memories_CR(hi), 0.1);
            % Ensure CRi is within the range [0, 1] and Fi is within the range [0,2]
            Fi = min(max(Fi, 0), 1);
            CRi = min(max(CRi, 0), 1);
        end

        function result = meanWL(obj, elements, weights)
            % Calculate the weighted Lehmer mean of elements.
            % Lehmer mean is calculated as the weighted sum of the squares
            % divided by the weighted sum of the elements.
        
            % Inputs:
            % elements: array of elements
            % weights: array of corresponding weights
        
            % Outputs:
            % result: weighted Lehmer mean

            numerator = sum((elements .^ 2) .* weights);
            denominator = sum(elements .* weights);
            if denominator == 0
                result = 0;
            else
                result = numerator / denominator;
            end
        end

        function result = meanWA(obj, elements, weights)
            % Calculate the weighted arithmetic mean of elements.
            % This is the standard weighted mean.
            result = sum(elements .* weights) / sum(weights);
        end

        function update_history(obj, S_F, S_CR, k)
            weights = ones(1, numel(S_F));
            if ~isempty(S_F)
                obj.memories_F(k) = obj.meanWL(S_F, weights);
            end
            if ~isempty(S_CR)
                obj.memories_CR(k) = obj.meanWA(S_CR, weights);
            end
        end

        function [final_root,final_score] = GA_evaluation(obj,verbose,visual_properties)
            rng(obj.seed);
            population = obj.generate_points(obj.population_size,obj.boundaries,obj.seed);
            
            fitness = zeros(1, obj.population_size);
            for i = 1:obj.population_size
                fitness(i) = obj.objective_function(population(i, :));
            end
            [best_fitness, best_idx] = min(fitness);
            best = population(best_idx, :);
            % subpopulation = zeros(obj.num_per_subpopulation,obj.dim,obj.population_size);
            % for i = 1:obj.population_size
            %     subpopulation(:, :, i) = obj.subpopulating(population(i, :), population, obj.num_per_subpopulation);
            % end
            
            % Animation Saving Setting
            if visual_properties.save_visual == true
                writerObj = VideoWriter(visual_properties.file_name);
                writerObj.FrameRate = 5;  % Adjust the frame rate as needed
                open(writerObj);
            end

            % Create a figure with visibility off
            if visual_properties.show_visual == false
                fig = figure('Visible', 'off');
            else 
                fig = figure('Visible', 'on');
                [X_surf, Y_surf] = meshgrid(linspace(obj.boundaries(1,1),obj.boundaries(1,2),100),linspace(obj.boundaries(2,1),obj.boundaries(2,2),100));
                XY_surf = [X_surf(:),Y_surf(:)];
                Z_surf = zeros(size(XY_surf,1),1);
            end

            memory_id=1;
            
            for gen = 1:obj.max_generation 
                if visual_properties.show_visual == true || visual_properties.save_visual == true
                    answ = [-6.437160, 0.155348;
                            -0.932122, 1.067870;
                            -0.155283, 6.439840;
                             0.163333, 6.122430;
                             0.667121, 0.690103;
                            -6.117110, -0.163476];
                    
                    for i=1:size(XY_surf,1)
                        Z_surf(i) = obj.repulsion_function(XY_surf(i,:));
                    end
                    Z_surf = reshape(Z_surf, size(X_surf));
                    subplot(2, 2, 1);
                    surf(X_surf, Y_surf, Z_surf);
                    view(0, 0);
                    title("Tampak Samping")

                    subplot(2, 2, 2);
                    surf(X_surf, Y_surf, Z_surf);
                    view(0, 90);
                    title("Tampak Atas")

                    subplot(2, 2, 3);
                    surf(X_surf, Y_surf, Z_surf);
                    view(0, -90);
                    title("Tampak Bawah")

                    subplot(2, 2, 4);
                    scatter(answ(:,1),answ(:,2), 50,'magenta',"LineWidth",2);
                    hold on;
                    scatter(population(:,1), population(:,2), 5, 'filled', 'blue');    
                    hold off;
                    rectangle('Position',[obj.boundaries(:,1)',(obj.boundaries(:,2)-obj.boundaries(:,1))'],'EdgeColor','#FF0000')
                    xlim(0.5*obj.boundaries(1,:));
                    ylim(0.5*obj.boundaries(2,:));
                end
                for ind=1:obj.population_size
                    obj.archive = obj.update_archive(population(ind,:),obj.archive);
                end
                S_F = [];
                S_CR = [];
                for i = 1:2:obj.population_size
                    [F_i, CR_i] = obj.update_parameter();
                    x_i = population(i, :);
                    parent = zeros(2, obj.dim); % because of one-point crossover
                
                    for j = 0:size(parent, 1)-1
                        subpopulation_i = obj.subpopulating(population(i+j, :), population, obj.num_per_subpopulation);
                        fitness_subpopulation_i = zeros(1,obj.num_per_subpopulation);
                        for k = 1:size(subpopulation_i,1)
                            fitness_subpopulation_i(k) = obj.objective_function(subpopulation_i(k,:));
                        end
                        [selected_subpopulation,selected_indices] = obj.selection(subpopulation_i,fitness_subpopulation_i);
                        [~,min_idx_parent] = min(fitness_subpopulation_i(selected_indices));
                        parent(j+1,:) = selected_subpopulation(min_idx_parent,:);
                    end
                    offspring_set = obj.crossover(parent(1,:),parent(2,:));
                   
                    for k=1:size(offspring_set,1)
                        trial = obj.mutate(offspring_set(k,:),F_i,obj.boundaries);
                        trial_fitness = obj.fitness_function(trial);
                        [id_closest_trial,closest_trial] = obj.closest_solution(trial,population);
                        closest_trial_fitness = obj.fitness_function(closest_trial);
                        if trial_fitness < closest_trial_fitness
                            fitness(id_closest_trial) = trial_fitness;
                            population(id_closest_trial,:) = trial;
                            S_F = [S_F, F_i];
                            S_CR = [S_CR, CR_i];
                            if trial_fitness < best_fitness
                                best_idx = i;
                                best = trial;
                                best_fitness = trial_fitness
                            end
                        end
                    end
                    if ~isempty(S_F) && ~isempty(S_CR)
                        obj.update_history(S_F, S_CR, memory_id);
                        memory_id = memory_id + 1;
                        if memory_id > obj.max_memories_size
                            memory_id = 1;
                        end
                    end
                end
                if verbose
                    fprintf("=========Generation %d=========\n", gen);
                    disp("Archive:");
                    disp(obj.archive);
                end

                if visual_properties.show_visual == true || visual_properties.save_visual == true
                    % Adjust aspect ratio
                    axis equal;
                    pbaspect([diff(xlim()) diff(ylim()) 1]);
        
                    % Maximize figure window
                    set(gcf, 'WindowState', 'maximized');
                    if ~isempty(obj.archive)
                        hold on
                        plot(obj.archive(:,1), obj.archive(:,2),'*',Color='magenta',MarkerSize=50);
                        hold off
                    end
                    title(sprintf("Generation %d",gen))
                    pause(0.05)
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(writerObj, frame);
                    end
                end
            end
            final_root = obj.archive;
            final_score = zeros(1, size(final_root,1));
            for fin_iter = 1:size(final_root,1)
                final_score(fin_iter) = obj.objective_function(final_root(fin_iter, :));
            end
        end

    end
end

