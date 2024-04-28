classdef SDGA < handle
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
        mutation_rate
        cluster_radius
        cluster_center
        cluster_iter_points
    end
    
    methods
        function obj = SDGA(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_rate,seed)
            obj.boundaries = boundaries;
            obj.dim = size(obj.boundaries,1);
            obj.m_cluster = m_cluster;
            obj.k_cluster = k_cluster;
            obj.m = m;
            obj.k_max = k_max;
            obj.epsilon = epsilon;
            obj.delta = delta;
            obj.gamma = gamma;
            obj.mutation_rate = mutation_rate;
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
                obj.cluster_radius = [obj.cluster_radius; norm(y - xt)];
            elseif (Fxt < Fy) && (Fxt < Fxc)
                obj.cluster_center = [obj.cluster_center; y];
                obj.cluster_radius = [obj.cluster_radius; norm(y - xt)];
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
                        plot(obj.cluster_center(c,1), obj.cluster_center(c,2), 'o','MarkerSize',1);
                        viscircles(obj.cluster_center(c,:), obj.cluster_radius(c), Color="magenta",LineWidth=3);
                        hold off
                    end
                    pause(0.25)
                    
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(writerObj, frame);
                    end
                end

                obj.cluster_iter_points = obj.recombination(obj.cluster_iter_points,obj.mutation_rate,obj.boundaries);
            
                k = k + 1;
            end

            if visual_properties.save_visual == true
                close(writerObj);
            end

        end

        function clean_roots = root_elimination(obj,root_archive)
            eligible_roots = [];
            for i = 1:size(root_archive,1)
                if obj.objective_function(root_archive(i,:)) < -1 + obj.epsilon    *1e4
                    eligible_roots = [eligible_roots; root_archive(i,:)];
                end
            end
            if size(eligible_roots,1) >1
                id_duplicated_roots = [];
                for i = 1:size(eligible_roots,1)
                    for j = i+1:size(eligible_roots,1)
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

        function [final_root,final_score] = GA_evaluation(obj,verbose, superverbose,visual_properties)
            obj.initialization();
            obj.clustering(visual_properties);
            if verbose == true
                fprintf("%d roots found!\n",numel(obj.cluster_radius))
            end
            
            archive = [];
            score = [];
            
            for i = 1:length(obj.cluster_center)
                subbound = [obj.cluster_center(i,:)' - obj.cluster_radius(i), obj.cluster_center(i,:)' + obj.cluster_radius(i)];
                [root, root_score] = obj.GA(obj.m,subbound,obj.k_max,obj.mutation_rate,obj.seed,superverbose);
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

