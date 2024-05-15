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
        cluster_iter_bestsol
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

        function [pop,BestSol] = generate_points(obj,npoint,boundaries,seed)
            rng(seed)
            dimension = size(boundaries,1);
            p = sobolset(dimension);
            p = scramble(p,'MatousekAffineOwen');
            A = net(p,npoint);
            empty_individual.Position = zeros(1,dimension);
            empty_individual.Cost = [];
            BestSol.Cost = inf;
            pop = repmat(empty_individual, npoint, 1);
            points = zeros(1,dimension);
            for i=1:npoint
                for j=1:dimension
                   points(:,j)=round((boundaries(j,1)+(boundaries(j,2)-boundaries(j,1)).*A(i,j))*100)/100;
                end
                pop(i).Position = points;
                pop(i).Cost = obj.objective_function(pop(i).Position);
                if pop(i).Cost<BestSol.Cost
                    BestSol = pop(i);
                end
            end
        end

        function [population,BestSol] = DE(obj,population,BestSol,boundaries,MaxIt,F,pCR,verbose)
            nPop = size(population,1);
            dimension = size(population(1).Position,2);
            VarSize = [1 dimension];
            BestCost = zeros(MaxIt, 1);
            for it = 1:MaxIt
                for i = 1:nPop

                    x = population(i).Position;

                    A = randperm(nPop);

                    A(A == i) = [];

                    a = A(1);
                    b = A(2);
                    c = A(3);

                    % Mutation
                    beta = repmat(F,VarSize);
                    dv_i = population(a).Position+beta.*(population(b).Position-population(c).Position);
                    dv_i = max(dv_i, boundaries(:,1)');
		            dv_i = min(dv_i, boundaries(:,2)');

                    % Crossover
                    trial = zeros(size(x));
                    j0 = randi([1 numel(x)]);
                    for j = 1:numel(x)
                        if j == j0 || rand <= pCR
                            trial(j) = dv_i(j);
                        else
                            trial(j) = x(j);
                        end
                    end

                    NewSol.Position = trial;
                    NewSol.Cost = obj.objective_function(NewSol.Position);

                    if NewSol.Cost<population(i).Cost
                        population(i) = NewSol;

                        if population(i).Cost<BestSol.Cost
                           BestSol = population(i);
                        end
                    end

                end

                % Update Best Cost
                BestCost(it) = BestSol.Cost;

                if verbose == true
                    % Show Iteration Information
                    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
                end

            end
        end


        function obj = initialization(obj)
            [population,bestsol] = obj.generate_points(obj.m_cluster,obj.boundaries,obj.seed);
            
            first_radius = (obj.boundaries(:, 2) - obj.boundaries(:, 1)) / 2;
            radii = min(first_radius);

            obj.cluster_center = bestsol.Position;
            obj.cluster_radius = radii;
            obj.cluster_iter_points = population;
            obj.cluster_iter_bestsol = bestsol;
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
                    func_point = obj.cluster_iter_points(i).Cost;
            
                    % If F(x_i)<gamma and x_i is not the center of existing cluster, x_i may have a possibility to become a cluster center
                    exist_in_cluster_center = any(vecnorm(obj.cluster_iter_points(i).Position - obj.cluster_center) < obj.epsilon);
                    if func_point < obj.gamma && exist_in_cluster_center == 0
                        potential_cluster_center = [potential_cluster_center; obj.cluster_iter_points(i).Position];
                    end
                end
                % Apply function cluster
                for i = 1:size(potential_cluster_center,1)
                    obj.function_cluster(potential_cluster_center(i,:));
                end
            
                if visual_properties.show_visual == true || visual_properties.save_visual == true
                    populationcell = arrayfun(@(x) x.Position, obj.cluster_iter_points, 'UniformOutput', false);
                    population = cell2mat(populationcell);
                    scatter(population(:,1), population(:,2), 30, 'filled');
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
                    pause(0.1)
                    
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(writerObj, frame);
                    end
                end

                % [obj.cluster_iter_points,~] = obj.reproduction(obj.cluster_iter_points, ...
                %     obj.boundaries,obj.mutation_factor,obj.crossover_rate,obj.seed);
                [obj.cluster_iter_points,obj.cluster_iter_bestsol] = obj.DE(obj.cluster_iter_points,obj.cluster_iter_bestsol,obj.boundaries,1,obj.mutation_factor,obj.crossover_rate,false);
            
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
                % [root, root_score] = obj.DE(subbound,obj.m,obj.k_max,obj.mutation_factor,obj.crossover_rate,obj.seed,superverbose);
                [popbound,BestSol] = obj.generate_points(obj.m,subbound,obj.seed);
                [~,BestSol] = obj.DE(popbound,BestSol,subbound,obj.k_max,obj.mutation_factor,obj.crossover_rate,superverbose);
                archive = [archive;BestSol.Position];
                score = [score;BestSol.Cost];
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

