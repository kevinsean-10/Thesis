classdef CDE < handle
    properties
        boundaries % boundaries of the NES
        dim % dimension of the NES
        seed % random state
        epsilon % accuracy
        criterion % criterion of the roots (based on the objective function used)
        gamma % ‘cut-off’ parameter for the objective function in clustering phase
        delta % minimum distance of solution
        m_cluster % number of individuals per generation in clustering phase
        k_cluster % maximum iteration in clustering phase
        m % number of individuals per generation in optimization phase
        k_max % maximum iteration in optimization phase
        mutation_factor % scale factor constant
        crossover_rate % crossover rate constant
        cluster_radius % set of cluster radius
        cluster_center % set of cluster radius
        cluster_iter_points % set of points in clustering phase
        cluster_iter_bestsol % set of best point in value in clustering phase
        fig % animation parameter
        writerObj % animation parameter
    end
    
    methods
        % INITIALIZATION
        function obj = CDE(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_factor,crossover_rate,seed)
            obj.boundaries = boundaries;
            obj.dim = size(obj.boundaries,1);
            obj.m_cluster = m_cluster;
            obj.k_cluster = k_cluster;
            obj.m = m;
            obj.k_max = k_max;
            obj.epsilon = epsilon;
            obj.criterion = epsilon;
            obj.delta = delta;
            obj.gamma = gamma;
            obj.mutation_factor = mutation_factor;
            obj.crossover_rate = crossover_rate;
            obj.seed = seed;
            obj = initialization(obj);
        end

        % % Problem 1
        % function F_array = system_equations(obj,x)
        %     f1 = exp(x(1)-x(2)) - sin(x(1)+x(2));
        %     f2 = (x(1)*x(2))^2 - cos(x(1)+x(2));
        %     F_array = [f1; f2];
        % end

        % % Problem 2
        % function F_array = system_equations(obj,x)
        %     f1 = 0.5 * sin(x(1) * x(2)) - 0.25 * x(2) / pi - 0.5 * x(1);
        %     f2 = (1 - 0.25 / pi) * (exp(2 * x(1)) - exp(1)) + exp(1) * x(2) / pi - 2 * exp(1) * x(1);
        %     F_array = [f1; f2];
        % end

        % Problem 7
        function F_array = system_equations(obj,x)
            f1 = x(1)^2-x(1)-x(2)^2-x(2)+x(3)^2;
            f2 = sin(x(2)-exp(x(1)));
            f3 = x(3)-log(abs(x(2)));
            F_array = [f1; f2; f3];
        end

        % Objective Function
        % INPUT
        % x: vectors with dim size
        % OUTPUT
        % res: real number
        function res = objective_function(obj, x)
            F_array = obj.system_equations(x);
            res = sum(abs(F_array));
            res = -1 / (1 + res);
            obj.criterion = -1 + obj.epsilon;
        end

        % Generate population
        % INPUT
        % npoint: desired population size
        % boundaries
        % seed
        % OUTPUT
        % pop: structured array of position and cost of the population
        % BestSol: structured array of position and cost of the best
        % solution
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

        % Mutation process of DE
        % INPUT
        % population: structured array
        % F: scale factor of mutation
        % VarSize: dimension of population
        % boundaries
        % permvec: randomized indices of population
        % OUTPUT
        % donor_vec: donor/mutate vector
        function donor_vec = mutate(obj,population,F,VarSize,boundaries, permvec)
            a = permvec(1);
            b = permvec(2);
            c = permvec(3);

            % Mutation
            beta = repmat(F,VarSize);
            donor_vec = population(permvec(1)).Position+beta.*(population(permvec(2)).Position-population(permvec(3)).Position);
            donor_vec = max(donor_vec, boundaries(:,1)');
            donor_vec = min(donor_vec, boundaries(:,2)');
        end

        % Crossover process of DE
        % INPUT
        % original_vec: original vector
        % donor_vec: donor vector
        % crossover_rate: crossover rate
        % OUTPUT
        % trial: trial vector
        function trial = crossover(obj,original_vec,donor_vec,crossover_rate)
            trial = zeros(size(original_vec));
            j0 = randi([1 numel(original_vec)]);
            for j = 1:numel(original_vec)
                if j == j0 || rand <= crossover_rate
                    trial(j) = donor_vec(j);
                else
                    trial(j) = original_vec(j);
                end
            end
        end

        % Evaluation using DE to get the minimum solution
        % INPUT
        % population: structured array
        % BestSol: structured array of previously found best position and cost
        % boundaries
        % MaxIt: maximum iteration
        % F: scale factor
        % pCR: crossover rate
        % verbose: print properties
        % visual_properties: visualization properties
        % OUTPUT
        % population: structured array
        % BestSol: structured array of best position and cost
        function [population,BestSol] = DE(obj,population,BestSol,boundaries,MaxIt,F,pCR,verbose,visual_properties)
            nPop = size(population,1);
            dimension = size(population(1).Position,2);
            VarSize = [1 dimension];
            BestCost = zeros(MaxIt, 1);
            for it = 1:MaxIt
                if visual_properties.show_visual == true || visual_properties.save_visual == true
                    pop_array = reshape([population.Position], dimension, [])';
                    hold on
                    pop_scatter = scatter(pop_array(:,1), pop_array(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.3,'MarkerEdgeAlpha',0,'MarkerFaceColor','blue');
                    hold off
                    pause(0.005)
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(obj.writerObj, frame);
                    end
                    delete (pop_scatter)
                end
                for i = 1:nPop

                    x = population(i).Position;

                    permvec = randperm(nPop); % Randomize the indices of the population
                    permvec(permvec == i) = [];

                    % Mutate
                    dv_i = obj.mutate(population,F,VarSize,boundaries,permvec);

                    % Crossover
                    trial = obj.crossover(x,dv_i,pCR);

                    % Offspring
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

        % Initialization of the variables
        function obj = initialization(obj)
            [population,bestsol] = obj.generate_points(obj.m_cluster,obj.boundaries,obj.seed);
            
            first_radius = (obj.boundaries(:, 2) - obj.boundaries(:, 1)) / 2;
            radii = min(first_radius);

            obj.cluster_center = bestsol.Position;
            obj.cluster_radius = radii;
            obj.cluster_iter_points = population;
            obj.cluster_iter_bestsol = bestsol;
        end

        % Function to generate new cluster or update the old one
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

        % Clustering phase
        % INPUT
        % visual_properties: visualization properties
        function clustering(obj,visual_properties)
            if visual_properties.save_visual == true
                obj.writerObj = VideoWriter(visual_properties.file_name);
                obj.writerObj.FrameRate = 5;  % Adjust the frame rate as needed
                open(obj.writerObj);
            end

            % Create a figure with visibility off
            if visual_properties.show_visual == false
                obj.fig = figure('Visible', 'off');
            else 
                obj.fig = figure('Visible', 'on');
                grid on;
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
                    population_array = reshape([obj.cluster_iter_points.Position], obj.dim, [])';
                    cluster_scatter = scatter(population_array(:,1), population_array(:,2), 36, 'filled', 'MarkerFaceColor','#F4A460', 'MarkerFaceAlpha', 0.3,'MarkerEdgeAlpha',0);
                    rectangle('Position', [obj.boundaries(1,1), obj.boundaries(2,1), obj.boundaries(1,2) - obj.boundaries(1,1), obj.boundaries(2,2) - obj.boundaries(2,1)], 'EdgeColor', '#FF0000', 'LineWidth', 5);

                    xlim(1.5 * obj.boundaries(1,:));
                    ylim(1.5 * obj.boundaries(2,:));
                    
                    % Adjust aspect ratio
                    axis equal;
                    pbaspect([diff(xlim()) diff(ylim()) 1]);
                    
                    % Maximize figure window
                    set(gcf, 'WindowState', 'maximized');
                    
                    for c = 1:size(obj.cluster_center,1)
                        hold on
                        % plot(obj.cluster_center(c,1), obj.cluster_center(c,2), 'o');
                        viscircles(obj.cluster_center(c,:), obj.cluster_radius(c), 'Color',"#008000", 'Linestyle','-','Linewidth',1);
                        hold off
                    end
                    pause(0.1)
                    
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(obj.writerObj, frame);
                    end
                end

                % [obj.cluster_iter_points,~] = obj.reproduction(obj.cluster_iter_points, ...
                %     obj.boundaries,obj.mutation_factor,obj.crossover_rate,obj.seed);
                vp = struct('show_visual',false, ...
                                            'save_visual', false, ...
                                            'file_name', 'sdde.avi');
                [obj.cluster_iter_points,obj.cluster_iter_bestsol] = obj.DE(obj.cluster_iter_points,obj.cluster_iter_bestsol,obj.boundaries,1,obj.mutation_factor,obj.crossover_rate,false,vp);
                
                k = k + 1;
            end
            if visual_properties.show_visual == true || visual_properties.save_visual == true
                delete(cluster_scatter);
            end
        end

        % Root elimination of the previously found solution to eliminate the double roots
        % INPUT
        % root_archive: set of found solutions
        % OUTPUT
        % clean_roots: set of "cleaned" solutions
        function clean_roots = root_elimination(obj,root_archive)
            eligible_roots = [];
            for i = 1:size(root_archive,1)
                if obj.objective_function(root_archive(i,:)) < -1 + obj.epsilon
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

        % Evaluation of DE in the clusters that potentially contain root
        % INPUT
        % verbose: print properties of DE_evaluation
        % superverbose: print properties of DE
        % visual_properties: visualization properties
        % OUTPUT
        % final_root: location of solution with minimum score
        % final_score: final_root's score
        % num_cluster: number of clusters that potentially contain root
        function [final_root,final_score,num_cluster] = DE_evaluation(obj,verbose, superverbose,visual_properties)
            obj.initialization();
            obj.clustering(visual_properties);
            num_cluster = numel(obj.cluster_radius);
            if verbose == true
                fprintf("%d roots found!\n",numel(obj.cluster_radius))
            end
            
            archive = [];
            score = [];

            for i = 1:length(obj.cluster_center)
                subbound = [obj.cluster_center(i,:)' - obj.cluster_radius(i), obj.cluster_center(i,:)' + obj.cluster_radius(i)];
                [popbound,BestSol] = obj.generate_points(obj.m,subbound,obj.seed);
                [~,BestSol] = obj.DE(popbound,BestSol,subbound,obj.k_max,obj.mutation_factor,obj.crossover_rate,superverbose,visual_properties);
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
            if visual_properties.show_visual == true
                figure(obj.fig)
                hold on
                for j = 1: size(final_root,1)
                    plot(final_root(j,1), final_root(j,2),'*','Color','magenta','MarkerSize',20);
                    pause(0.25)
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(obj.writerObj, frame);
                    end
                end
                hold off
                if visual_properties.save_visual == true
                close(obj.writerObj);
                end
            end
        end

    end
end

