classdef HGA < handle
    properties
        boundaries
        dim
        population_size
        max_generation
        mutation_rate
        delta
        epsilon
        parts
        seed
        cluster
        archive
        score
        final_root
        final_score
    end
    methods
        function obj = HGA(boundaries, population_size, parts, max_generation,mutation_rate,epsilon, delta, seed)
            obj.boundaries = boundaries;
            obj.dim = size(boundaries,1);
            obj.population_size = population_size;
            obj.max_generation = max_generation;
            obj.parts = parts;
            obj.mutation_rate = mutation_rate;
            obj.epsilon = epsilon;
            obj.delta = delta;
            obj.seed = seed;
            obj.cluster = {};
            obj.archive = {};
            obj.score = [];
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

        function selected_population = selection(obj, population, fitness)
            pop_size = size(population, 1);
            selection_probs = 1 ./ (fitness + 1); % add one to avoid division by zero
            total_probs = sum(selection_probs);
            selection_probs = selection_probs / total_probs;
            selected_indices = randsample(1:pop_size, pop_size, true, selection_probs);
            selected_population = population(selected_indices, :);
        end

        function [offspring1, offspring2] = crossover(obj, parent1, parent2)
            dimension = length(parent1);
            crossover_point = randi([1, dimension], 1);
            offspring1 = [parent1(1:crossover_point), parent2(crossover_point+1:end)];
            offspring2 = [parent2(1:crossover_point), parent1(crossover_point+1:end)];
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

        function [best_individual, best_fitness] = GA(obj, population_size, boundaries, max_generation, mutation_rate, seed, print_stat)
            rng(seed);
            dimension = size(boundaries, 2);
            population = obj.generate_points(population_size, boundaries, seed);
            fitness = zeros(1, population_size);
            for i = 1:population_size
                fitness(i) = obj.objective_function(population(i, :));
            end
            [best_fitness, best_idx] = min(fitness);
            best_individual = population(best_idx, :);
        
            for generation = 1:max_generation
                selected_population = obj.selection(population, fitness);
                offspring_population = obj.recombination(selected_population, mutation_rate, boundaries);
                population = offspring_population;
                fitness = zeros(1, population_size);
                for i = 1:population_size
                    fitness(i) = obj.objective_function(population(i, :));
                end
                [best_fitness, best_idx] = min(fitness);
                best_individual = population(best_idx, :);
                if print_stat == true
                    disp("=========Generation " + generation + "=========");
                    disp("Best Individual: ");
                    disp(best_individual);
                    disp("Best Score: ");
                    disp(best_fitness);
                    disp(" ");
                end
            end
        end

        function res = slice_hypercube(obj, lower_bounds, upper_bounds, parts)
            interval = (upper_bounds-lower_bounds)/(parts);
            dimension = numel(lower_bounds);
            
            % Create a cell array of arrays, each containing points spaced interval apart for each dimension
            points = cell(1, dimension);
            for i = 1:dimension
                points{i} = linspace(lower_bounds(i), upper_bounds(i)-interval(i), parts);
            end
            
            % Use ndgrid to create a grid of points in n-dimensional space
            [grids{1:dimension}] = ndgrid(points{:});
            
            % Combine the grid points into a single matrix
            grid_points = reshape(cat(dimension+1, grids{:}), [], dimension);
            
            % % Generate all vertices for smaller hypercubes within each grid cell
            zero_interval = cell(dimension,1);
            for i = 1:dimension
                zero_interval{i} = [0,interval(i)];
            end
            
            [offsets_grid{1:dimension}] = ndgrid(zero_interval{:});
            offsets_grid_points = reshape(cat(dimension+1, offsets_grid{:}), [], dimension);
            res = zeros([size(offsets_grid_points, 1),size(grid_points)]);
            for i = 1:size(offsets_grid_points, 1)
                res(i,:,:) = grid_points + offsets_grid_points(i,:);
            end
        end

        function clustering(obj)
            lower_bounds = obj.boundaries(:,1); upper_bounds = obj.boundaries(:,2); 
            interval = (upper_bounds-lower_bounds)/(obj.parts);
            
            hypercubes_edges = obj.slice_hypercube(lower_bounds,upper_bounds,obj.parts);
            
            for hypercube_id = 1:size(hypercubes_edges, 2)
                X0 = hypercubes_edges(:,hypercube_id,:);
                X0 = reshape(X0, [], size(X0, 3));
                
                F_list = zeros(obj.dim,size(X0,1));
                for i = 1:size(X0,1)
                    F_list(:,i) = obj.system_equations(X0(i,:));
                end
                
                product_combination = zeros(size(F_list, 1), nchoosek(size(F_list, 2), 2));
                for i = 1:size(F_list, 1)
                    combinations = nchoosek(F_list(i,:), 2);
                    product_combination(i,:) = prod(combinations, 2)';
                end
                
                change_sign = any(product_combination < 0, 2);
                if all(change_sign)
                    obj.cluster{end+1} = X0;
                end
            end
            % Convert cluster to array
            obj.cluster = cat(3, obj.cluster{:});
        end

        function final_root = root_elimination(obj,root_archive)
            eligible_roots = [];
            for i = 1:size(root_archive,1)
                if obj.objective_function(root_archive(i,:)) < -1 + obj.epsilon
                    eligible_roots = [eligible_roots; root_archive(i,:)];
                end
            end
            
            id_duplicated_roots = [];
            for i = 1:length(eligible_roots)
                for j = i+1:length(eligible_roots)
                    if norm(eligible_roots(i,:) - eligible_roots(j,:)) < obj.delta
                        id_duplicated_roots = [id_duplicated_roots; [i, j]];
                    end
                end
            end
            
            id_duplicated_roots = unique(id_duplicated_roots, 'rows');
            
            deselected_id_duplicated_roots = [];
            for i = 1:size(id_duplicated_roots, 1)
                root_a = obj.objective_function(eligible_roots(id_duplicated_roots(i, 1),:));
                root_b = obj.objective_function(eligible_roots(id_duplicated_roots(i, 2),:));
                if root_a <= root_b
                    id_duplicated_root = id_duplicated_roots(i, 2);
                else
                    id_duplicated_root = id_duplicated_roots(i, 1);
                end
                deselected_id_duplicated_roots = [deselected_id_duplicated_roots; id_duplicated_root];
            end
            
            if ~isempty(deselected_id_duplicated_roots)
                unique_roots = true(size(eligible_roots,1),1);
                unique_roots(deselected_id_duplicated_roots) = false;
                final_root = eligible_roots(unique_roots,:);
            else
                final_root = eligible_roots;
            end
        end
        
        function [final_root,final_score] = GA_evaluation(obj, verbose, superverbose)
            obj.clustering()
            if verbose == true
                fprintf("Number of clusters containing root: %d\n", size(obj.cluster, 3));
            end
            for i = 1:size(obj.cluster, 3)
                subbound = zeros(obj.dim, 2);
                for d = 1:size(obj.cluster, 2)
                    subbound(d, :) = [min(obj.cluster(:, d, i)), max(obj.cluster(:, d, i))];
                end
        
                [root, root_score] = obj.GA(obj.population_size,subbound,obj.max_generation,obj.mutation_rate,obj.seed,superverbose);
                
                obj.archive{end+1} = root;
                obj.score(end+1) = root_score;
                
                if verbose
                    fprintf('\n====== Cluster %d ======\n', i);
                    disp('Roots =');
                    disp(obj.archive);
                end
            end

            % Convert cell array to matrix
            archive_matrix = zeros(numel(obj.archive), numel(obj.archive{1}));
            for i = 1:numel(obj.archive)
                archive_matrix(i, :) = obj.archive{i};
            end

            final_root = root_elimination(obj, archive_matrix);
            final_score = zeros(1, size(final_root,1));
            for fin_iter = 1:size(final_root,1)
                final_score(fin_iter) = obj.objective_function(final_root(fin_iter, :));
            end
            obj.final_root = final_root;
            obj.final_score = final_score;
        end
        
        function visualization2D(obj,visual_properties)
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
            end

            % Adjust aspect ratio
            axis equal;
            pbaspect([diff(xlim()) diff(ylim()) 1]);
            
            % Maximize figure window
            set(gcf, 'WindowState', 'maximized');

            xlim(obj.boundaries(1,:));
            ylim(obj.boundaries(2,:));
            hold on
            rectangle('Position',[obj.boundaries(:,1)',(obj.boundaries(:,2)-obj.boundaries(:,1))'],'EdgeColor','#FF0000')
            for i = 1:size(obj.cluster, 3)
                subbound = zeros(obj.dim, 2);
                for d = 1:size(obj.cluster, 2)
                    subbound(d, :) = [min(obj.cluster(:, d, i)), max(obj.cluster(:, d, i))];
                end
                
                rectangle('Position',[subbound(:,1)',(subbound(:,2)-subbound(:,1))'],'EdgeColor','#4DBEEE')
                pause(0.25)
                if visual_properties.save_visual == true
                    frame = getframe(gcf);
                    writeVideo(writerObj, frame);
                end
            end

            for j = 1: size(obj.final_root,1)
                plot(obj.final_root(j,1), obj.final_root(j,2), '*',Color='magenta',MarkerSize=50);
                pause(0.25)
                if visual_properties.save_visual == true
                    frame = getframe(gcf);
                    writeVideo(writerObj, frame);
                end
            end
            hold off

            if visual_properties.save_visual == true
                close(writerObj);
            end

        end


    end
end