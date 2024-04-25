classdef SDDE < handle
    properties
        boundaries
        dim
        seed
        epsilon
        gamma
        delta
        m_cluster
        k_cluster
        m
        k_max
        mutation_factor
        crossover_rate
        cluster_radius
        cluster_center
        cluster_iter_points
    end
    
    methods
        function obj = SDDE(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_factor,crossover_rate,seed)
            obj.boundaries = boundaries;
            obj.dim = size(obj.boundaries,1);
            obj.m_cluster = m_cluster;
            obj.k_cluster = k_cluster;
            obj.m = m;
            obj.k_max = k_max;
            obj.epsilon = epsilon;
            obj.delta = delta;
            obj.gamma = gamma;
            obj.mutation_factor = mutation_factor;
            obj.crossover_rate = crossover_rate;
            obj.seed = seed;
            obj = initialization(obj);
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
        
            % Generate three distinct individuals xr1, xr1, xr1 from the current population randomly
            pop_ids = 1:size(population, 1);
            if nargin > 4 && ~isempty(x_i_id)
                index_to_delete = x_i_id;
            else
                index_to_delete = find(all(population_copy == x_i, 2)); % Ensure that x_i is excluded from the selected subpopulation
            end
            pop_ids_no_i = setdiff(pop_ids, index_to_delete);
            population = population(pop_ids_no_i, :);
        
            % Mutation form the donor/mutation vector
            dv_i = obj.mutate(population, mutation_factor);
        
            % Set penalty for every donor vector that violates the boundaries
            % for j = 1:size(dv_i, 2)
            %     if dv_i(j) < boundaries(j, 1)
            %         dv_i(j) = (x_i(j) + boundaries(j, 1)) / 2;
            %     elseif dv_i(j) > boundaries(j, 2)
            %         dv_i(j) = (x_i(j) + boundaries(j, 2)) / 2;
            %     end
            % end
        end

        function result_population = reproduction(obj,population,boundaries,mutation_factor,crossover_rate,seed)
            rng(seed); % Set random seed
            for i = 1:size(population,1)
                x_i = population(i,:);
                dv_i = obj.mutation_penalty(x_i, population, boundaries, mutation_factor,i);
                trial = obj.crossover(x_i, dv_i, crossover_rate);
                if obj.objective_function(trial)<=obj.objective_function(x_i)
                    population(i,:) = trial;
                end
            end
            fitness_pop = zeros(1,size(population,1));
            for i = 1:size(population,1)
                fitness_pop(i) = obj.objective_function(population(i,:));
            end
            [~,id_fit] = sort(fitness_pop);
            result_population = population(id_fit,:);
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

        function obj = initialization(obj)
            points = obj.generate_points(obj.m_cluster,obj.boundaries,obj.seed);
            fitness = zeros(size(points,1),1);
            for i=1:size(points,1)
                fitness(i)=obj.objective_function(points(i,:));
            end
            [~, best_idx] = min(fitness);
            x_prime = points(best_idx, :);
            
            first_radius = (obj.boundaries(:, 2) - obj.boundaries(:, 1)) / 2;
            radii = min(first_radius);

            obj.cluster_center = x_prime;
            obj.cluster_radius = radii;
            obj.cluster_iter_points = points;
        end

        function function_cluster(obj,y)
            dist_list = vecnorm(obj.cluster_center - y, 2, 2);
            [~, min_dist_id] = min(dist_list);
            xc = obj.cluster_center(min_dist_id, :);
            xt = (xc + y) / 2;

            Fxt = obj.objective_function(xt');
            Fxc = obj.objective_function(xc');
            Fy = obj.objective_function(y');

            if (Fxt > Fy) && (Fxt > Fxc)
                obj.cluster_center = [obj.cluster_center; y];
                obj.cluster_radius = [obj.cluster_radius; norm((y - xt),2)];
            elseif (Fxt < Fy) && (Fxt < Fxc)
                obj.cluster_center = [obj.cluster_center; y];
                obj.cluster_radius = [obj.cluster_radius; norm((y - xt),2)];
                obj.function_cluster(xt);
            elseif Fy < Fxc
                obj.cluster_center(min_dist_id, :) = y;
            end

            obj.cluster_radius(min_dist_id) = norm(y - xt);
        end

        function clustering(obj,visual_properties)
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

            k = 0;
            while k < obj.k_cluster
                potential_cluster_center = [];
                for i = 1:size(obj.cluster_iter_points,1)
                    func_point = obj.objective_function(obj.cluster_iter_points(i,:));
            
                    % If F(x_i)<gamma and x_i is not the center of existing cluster, x_i may have a possibility to become a cluster center
                    exist_in_cluster_center = any(vecnorm(obj.cluster_iter_points(i,:) - obj.cluster_center) < obj.epsilon);
                    if func_point < obj.gamma && exist_in_cluster_center == 0
                        potential_cluster_center = [potential_cluster_center; obj.cluster_iter_points(i,:)];
                    end
                end
                % Apply function cluster
                for i = 1:size(potential_cluster_center,1)
                    obj.function_cluster(potential_cluster_center(i,:));
                end
            
                if visual_properties.show_visual == true || visual_properties.save_visual == true
                    scatter(obj.cluster_iter_points(:,1), obj.cluster_iter_points(:,2), 30, 'filled');
                    rectangle('Position',[obj.boundaries(:,1)',(obj.boundaries(:,2)-obj.boundaries(:,1))'],'EdgeColor','#FF0000')
                    xlim(1.5 * obj.boundaries(1,:));
                    ylim(1.5 * obj.boundaries(2,:));
                    
                    % Adjust aspect ratio
                    axis equal;
                    pbaspect([diff(xlim()) diff(ylim()) 1]);
                    
                    % Maximize figure window
                    set(gcf, 'WindowState', 'maximized');
                    
                    for c = 1:size(obj.cluster_center,1)
                        hold on
                        plot(obj.cluster_center(c,1), obj.cluster_center(c,2), 'o');
                        viscircles(obj.cluster_center(c,:), obj.cluster_radius(c), Color="#00FFFF", Linestyle='-.');
                        hold off
                    end
                    pause(0.25)
                    
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(writerObj, frame);
                    end
                end

                obj.cluster_iter_points = obj.reproduction(obj.cluster_iter_points, ...
                    obj.boundaries,obj.mutation_factor,obj.crossover_rate,obj.seed);
            
                k = k + 1;
            end

            if visual_properties.save_visual == true
                close(writerObj);
            end

        end

        function clean_roots = root_elimination(obj,root_archive)
            eligible_roots = [];
            for i = 1:size(root_archive,1)
                if obj.objective_function(root_archive(i,:)) < -1 + obj.epsilon
                    eligible_roots = [eligible_roots; root_archive(i,:)];
                end
            end
            
            if size(eligible_roots,1) >1
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
                    clean_roots = eligible_roots(unique_roots,:);
                else
                    clean_roots = eligible_roots;
                end
            else
                clean_roots = eligible_roots;
            end
        end

        function [final_root,final_score] = DE_evaluation(obj,verbose, superverbose,visual_properties)
            obj.initialization();
            obj.clustering(visual_properties);
            if verbose == true
                fprintf("%d roots found!\n",numel(obj.cluster_radius))
            end
            
            archive = [];
            score = [];
            
            for i = 1:length(obj.cluster_center)
                subbound = [obj.cluster_center(i,:)' - obj.cluster_radius(i), obj.cluster_center(i,:)' + obj.cluster_radius(i)];
                [root, root_score] = obj.DE(subbound,obj.m,obj.k_max,obj.mutation_factor,obj.crossover_rate,obj.seed,superverbose);
                archive = [archive;root];
                score = [score;root_score];
                if verbose == true
                    fprintf('\n====== Cluster %d ======\n', i);
                    disp(archive);
                end
            end
            final_root = obj.root_elimination(archive);
            final_score = zeros(1, size(final_root,1));
            for fin_iter = 1:size(final_root,1)
                final_score(fin_iter) = obj.objective_function(final_root(fin_iter, :));
            end
        end


    end
end

