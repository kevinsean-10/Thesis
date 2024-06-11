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

        % % Problem 4
        % function F_array = system_equations(obj,x)
        %     % g1 = x(1)*x(2)^3 / 12 - (x(1) - 2*x(3))*(x(2) - 2*x(3))^3 / 12 - 9369;
        %     % g2 = 2*(x(2) - x(3))^2 * (x(1) - x(3))^2 * x(3) / (x(2) + x(1) - 2*x(3)) - 6835;
        %     f1 = x(1)*x(2) - (x(1) - 2*x(3))*(x(2) - 2*x(3)) - 165;
        %     f2 = x(1)*x(2)^3 / 12 - (x(1) - 2*x(3))*(x(2) - 2*x(3))^3 / 12 - 9369;
        %     f3 = (2*((x(2)-x(3))^2)*((x(1)-x(3))^2)*x(3))/(x(2)+x(1)-2*x(3)) - 6835;
        %     F_array = [f1; f2; f3];
        % end

        % % Problem 5
        % function F_array = system_equations(obj,x)
        %     f1 = 2*x(1) + x(2) + x(3) + x(4) + x(5) - 6;
        %     f2 = x(1) + 2*x(2) + x(3) + x(4) + x(5) - 6;
        %     f3 = x(1) + x(2) + 2*x(3) + x(4) + x(5) - 6;
        %     f4 = x(1) + x(2) + x(3) + 2*x(4) + x(5) - 6;
        %     f5 = x(1)*x(2)*x(3)*x(4)*x(5) - 1;
        %     F_array = [f1; f2;f3;f4;f5];
        % end
        
        % Problem 7
        function F_array = system_equations(obj,x)
            f1 = x(1)^2-x(1)-x(2)^2-x(2)+x(3)^2;
            f2 = sin(x(2)-exp(x(1)));
            f3 = x(3)-log(abs(x(2)));
            F_array = [f1; f2; f3];
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

        function donor_vec = mutate(obj,population,F,VarSize,boundaries, permvec)
            a = permvec(1);
            b = permvec(2);
            c = permvec(3);

            % Mutation
            beta = repmat(F,VarSize);
            donor_vec = population(a).Position+beta.*(population(b).Position-population(c).Position);
            donor_vec = max(donor_vec, boundaries(:,1)');
            donor_vec = min(donor_vec, boundaries(:,2)');
        end

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

        function [population,BestSol] = DE(obj,population,BestSol,boundaries,MaxIt,F,pCR,verbose)
            nPop = size(population,1);
            dimension = size(population(1).Position,2);
            VarSize = [1 dimension];
            BestCost = zeros(MaxIt, 1);
            for it = 1:MaxIt
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

        function clean_roots = root_elimination(obj,root_archive)
            eligible_roots = [];
            for i = 1:size(root_archive,1)
                if obj.objective_function(root_archive(i,:)) < -1 + obj.epsilon
                    eligible_roots = [eligible_roots; root_archive(i,:)];
                end
            end
            eligible_roots
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

        function [final_root,final_score] = DE_evaluation(obj, verbose, superverbose)
            obj.clustering()
            num_cluster = size(obj.cluster,1);
            if verbose == true
                fprintf("Number of clusters containing root: %d\n", size(obj.cluster, 3));
            end
            for i = 1:size(obj.cluster, 3)
                subbound = zeros(obj.dim, 2);
                for d = 1:size(obj.cluster, 2)
                    subbound(d, :) = [min(obj.cluster(:, d, i)), max(obj.cluster(:, d, i))];
                end
                [popbound,BestSol] = obj.generate_points(obj.population_size,subbound,obj.seed);
                [~,BestSol] = obj.DE(popbound,BestSol,subbound,obj.max_generation,obj.mutation_factor,obj.crossover_rate,superverbose);
                
                obj.archive{end+1} = BestSol.Position;
                obj.score(end+1) = BestSol.Cost;
                
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

            final_root = obj.root_elimination(archive_matrix);
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
                plot(obj.final_root(j,1), obj.final_root(j,2), '*','Color','magenta');
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

