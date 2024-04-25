classdef HDE < handle   
    properties
        boundaries
        dim
        population_size
        max_generation
        mutation_factor
        crossover_rate
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
        function obj = HDE(boundaries, population_size, parts, max_generation,mutation_factor, crossover_rate,epsilon, delta, seed)
            obj.boundaries = boundaries;
            obj.dim = size(boundaries,1);
            obj.population_size = population_size;
            obj.max_generation = max_generation;
            obj.parts = parts;
            obj.mutation_factor = mutation_factor;
            obj.crossover_rate = crossover_rate;
            obj.epsilon = epsilon;
            obj.delta = delta;
            obj.seed = seed;
            obj.cluster = {};
            obj.archive = {};
            obj.score = [];
            obj.final_root = [];
            obj.final_score = [];
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

        function mutant = mutate(obj, population, F)
            % Mutation function for DE
            % Vectorized mutation operation
            [~, indices] = sort(randperm(size(population, 1)));
            r = population(indices(1:3), :);
            mutant = r(1, :) + F * (r(2, :) - r(3, :));
        end

        function trial = crossover(obj, target, mutant, CR)
            % Crossover function for DE
            cross_points = rand(size(target)) < CR;
            % Ensure at least one true crossover point
            if ~any(cross_points(:))
                cross_points(randi(size(target, 2))) = true;
            end
            trial = mutant;
            trial(cross_points) = target(cross_points);
        end

        function dv_i = mutation_penalty(obj, x_i, population, boundaries, mutation_factor,x_i_id)
            % Mutation function with penalty for DE
            % Inputs:
            % x_i: target x_i
            % subpop_i: number of individuals closest to x_i
            % boundaries: boundaries/constraints of the function
            % scaling_factor: scaling factor of the function
            % Output:
            % dv_i: donor vector that has been mutated and penalized.
        
            % Generate three distinct individuals xr1, xr1, xr1 from 
            % the current population randomly
            population_copy = population;
            pop_ids = 1:size(population_copy, 1);
            if nargin > 4 && ~isempty(x_i_id)
                index_to_delete = x_i_id;
            else
                index_to_delete = find(all(population_copy == x_i, 2)); 
                % Ensure that x_i is excluded from the selected subpopulation
            end
            pop_ids_no_i = setdiff(pop_ids, index_to_delete);
            population_copy = population_copy(pop_ids_no_i, :);
        
            % Mutation form the donor/mutation vector
            dv_i = obj.mutate(population_copy, mutation_factor);
        
            % Set penalty for every donor vector that violates the boundaries
            for j = 1:size(dv_i, 2)
                if dv_i(j) < boundaries(j, 1)
                    dv_i(j) = (x_i(j) + boundaries(j, 1)) / 2;
                elseif dv_i(j) > boundaries(j, 2)
                    dv_i(j) = (x_i(j) + boundaries(j, 2)) / 2;
                end
            end
        end

        function [best_point, best_score] = DE(obj, boundaries,population_size,max_generation,mutation_factor,crossover_rate,seed,print_stat)
            rng(seed); % Set random seed
            population = obj.generate_points(population_size, boundaries, seed);
            fitness = zeros(1, population_size);
            for i = 1:population_size
                fitness(i) = obj.objective_function(population(i, :));
            end
            [best_score, best_idx] = min(fitness);
            best_point = population(best_idx,:);

            for gen = 1:max_generation
                for i = 1:population_size
                    x_i = population(i,:);
                    dv_i = obj.mutation_penalty(x_i, population, boundaries, mutation_factor,i);
                    trial = obj.crossover(x_i, dv_i, crossover_rate);
                    trial_fitness = obj.objective_function(trial);

                    if trial_fitness <= fitness(i)
                        fitness(i) = trial_fitness;
                        population(i,:) = trial;
                        if trial_fitness < fitness(best_idx)
                            best_idx = i;
                            best_point = trial;
                            best_score = fitness(best_idx);
                        end
                    end
                end

                if print_stat
                    fprintf("=========Generation %d=========\n", gen);
                    fprintf("Best Point: %s with score %.4f\n", mat2str(best_point), best_score);
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

        function [final_root,final_score] = DE_evaluation(obj, verbose, superverbose)
            obj.clustering()
            if verbose == true
                fprintf("Number of clusters containing root: %d\n", size(obj.cluster, 3));
            end
            for i = 1:size(obj.cluster, 3)
                subbound = zeros(obj.dim, 2);
                for d = 1:size(obj.cluster, 2)
                    subbound(d, :) = [min(obj.cluster(:, d, i)), max(obj.cluster(:, d, i))];
                end
        
                [root, root_score] = obj.DE(subbound,obj.population_size,obj.max_generation,obj.mutation_factor,obj.crossover_rate,obj.seed,superverbose);
                
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
                plot(obj.final_root(j,1), obj.final_root(j,2), '*',Color='magenta');
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

