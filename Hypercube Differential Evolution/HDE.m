classdef HDE < handle   
    properties
        boundaries % boundaries of the NES
        dim % dimension of the NES
        population_size % number of individuals per generation
        max_generation % maximum generation/iteration
        mutation_factor % scale factor constant
        crossover_rate % crossover rate constant
        delta % minimum distance of solution
        epsilon % accuracy
        parts % number of desired hypercubes in the boundaries (N)
        seed % random state
        cluster % set of choosen hypercubes
        archive % set of solutions
        score % set of solutions score
        final_root % set of "cleaned" solutions
        final_score % set of "cleaned" solutions score
        hypercubes_edges % set of all hypercubes
        fig % animation parameter
        writerObj % animation parameter
    end
    
    methods
        % INITIALIZATION
        function obj = HDE(boundaries, population_size, parts, ...
                max_generation,mutation_factor, crossover_rate,epsilon, ...
                delta, seed)
            obj.boundaries = boundaries;
            obj.dim = size(boundaries,1);
            obj.population_size = population_size;
            obj.max_generation = max_generation;
            obj.parts = parts;
            obj.mutation_factor = mutation_factor;
            obj.crossover_rate = crossover_rate;
            obj.hypercubes_edges = [];
            obj.epsilon = epsilon;
            obj.delta = delta;
            obj.seed = seed;
            obj.cluster = {};
            obj.archive = {};
            obj.score = [];
            obj.final_root = [];
            obj.final_score = [];
        end
        
        % Problem 1
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
        
        % % Problem 3
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
            donor_vec = population(a).Position+beta.*(population(b).Position-population(c).Position);
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

        % slicing the domain into hypercubes
        % INPUT
        % lower_bounds, upper_bounds: array of lower and upper bounds
        % parts: squared number to indicate the number of desired hypercubes in the domain
        % OUTPUT
        % res: set of hypercubes with size number of corner point of
        % a single hypercube by number of hypercubes by dimension
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

        % Selection process of hypercubes to identify which hypercube contains root
        function clustering(obj)
            lower_bounds = obj.boundaries(:,1); upper_bounds = obj.boundaries(:,2); 
            interval = (upper_bounds-lower_bounds)/(obj.parts);
            
            obj.hypercubes_edges = obj.slice_hypercube(lower_bounds,upper_bounds,obj.parts);
            
            for hypercube_id = 1:size(obj.hypercubes_edges, 2)
                X0 = obj.hypercubes_edges(:,hypercube_id,:);
                X0 = reshape(X0, [], size(X0, 3));
                
                F_list = zeros(obj.dim,size(X0,1));
                for i = 1:size(X0,1)
                    F_list(:,i) = obj.system_equations(X0(i,:));
                end
                
                % Selection process based of whether the f in NES changes sign or not
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

        % Root elimination of the previously found solution to eliminate the double roots
        % INPUT
        % root_archive: set of found solutions
        % OUTPUT
        % clean_roots: set of "cleaned" solutions
        function clean_roots = root_elimination(obj,root_archive)
            eligible_roots = [];
            % root criterion
            for i = 1:size(root_archive,1)
                if obj.objective_function(root_archive(i,:)) < -1 + obj.epsilon
                    eligible_roots = [eligible_roots; root_archive(i,:)];
                end
            end

            % eliminate double roots
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

        % Evaluation of DE in the hypercubes that potentially contain root
        % INPUT
        % verbose: print properties of DE_evaluation
        % superverbose: print properties of DE
        % visual_properties: visualization properties
        % OUTPUT
        % final_root: location of solution with minimum score
        % final_score: final_root's score
        % num_cluster: number of hypercubes that potentially contain root
        function [final_root,final_score,num_cluster] = DE_evaluation(obj, verbose, superverbose,visual_properties)
            if visual_properties.save_visual == true
                obj.writerObj = VideoWriter(visual_properties.file_name);
                obj.writerObj.FrameRate = 5;  % Adjust the frame rate as needed
                open(obj.writerObj);
            end

            % Create a figure with visibility off
            if visual_properties.show_visual == false
                fig = figure('Visible', 'off');
            else 
                fig = figure('Visible', 'on');
            end
            
            obj.clustering()
            num_cluster = size(obj.cluster,3);

            if visual_properties.show_visual == true || visual_properties.save_visual == true
                % axis equal;
                pbaspect([diff(xlim()) diff(ylim()) 1]);
                
                % Maximize figure window
                set(gcf, 'WindowState', 'maximized');
    
                xlim(obj.boundaries(1,:));
                ylim(obj.boundaries(2,:));
                title(sprintf('N=%d',obj.parts^2),"FontSize",24)
                hold on
                rectangle('Position',[obj.boundaries(:,1)',(obj.boundaries(:,2)-obj.boundaries(:,1))'],'EdgeColor','#FF0000','LineWidth',5)
                for i = 1:size(obj.hypercubes_edges, 2)
                    subbound = zeros(obj.dim, 2);
                    for d = 1:size(obj.hypercubes_edges, 3)
                        subbound(d, :) = [min(obj.hypercubes_edges(:, i,d)), max(obj.hypercubes_edges(:, i,d))];
                    end
                    rectangle('Position',[subbound(:,1)',(subbound(:,2)-subbound(:,1))'],'EdgeColor','#F4A460')
                end
                hold off
            end


            if verbose == true
                fprintf("Number of clusters containing root: %d\n", size(obj.cluster, 3));
            end
            for i = 1:size(obj.cluster, 3)
                subbound = zeros(obj.dim, 2);
                for d = 1:size(obj.cluster, 2)
                    subbound(d, :) = [min(obj.cluster(:, d, i)), max(obj.cluster(:, d, i))];
                end
                if visual_properties.show_visual == true || visual_properties.save_visual == true
                    rectangle('Position',[subbound(:,1)',(subbound(:,2)-subbound(:,1))'],'EdgeColor','#4DBEEE','LineWidth',3)
                    pause(0.25)
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(obj.writerObj, frame);
                    end
                end
                [popbound,BestSol] = obj.generate_points(obj.population_size,subbound,obj.seed);
                [~,BestSol] = obj.DE(popbound,BestSol,subbound,obj.max_generation,obj.mutation_factor,obj.crossover_rate,superverbose,visual_properties);
                
                obj.archive{end+1} = BestSol.Position;
                obj.score(end+1) = BestSol.Cost;

                % Convert cell array to matrix
                archive_matrix = reshape([obj.archive{:}], obj.dim, [])';
                if visual_properties.show_visual == true || visual_properties.save_visual == true
                    hold on
                    archive_scatter = scatter(archive_matrix(:,1),archive_matrix(:,2),'filled','MarkerEdgeAlpha',0,'MarkerFaceColor',[34, 65, 65]/255);
                    hold off
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(obj.writerObj, frame);
                    end
                end
                
                if verbose
                    fprintf('\n====== Cluster %d ======\n', i);
                    disp('Roots =');
                    disp(obj.archive);
                end
            end

            final_root = obj.root_elimination(archive_matrix);
            final_score = zeros(1, size(final_root,1));
            for fin_iter = 1:size(final_root,1)
                final_score(fin_iter) = obj.objective_function(final_root(fin_iter, :));
            end
            obj.final_root = final_root;
            obj.final_score = final_score;

            if visual_properties.show_visual == true || visual_properties.save_visual == true
                delete (archive_scatter);
                hold on
                plot(obj.final_root(:,1), obj.final_root(:,2),'*','Color','magenta','MarkerSize',20);      
                hold off
                pause(0.25)
                if visual_properties.save_visual == true
                    frame = getframe(gcf);
                    writeVideo(obj.writerObj, frame);
                end
            end
            if visual_properties.save_visual == true
                close(obj.writerObj);
            end
        end
    end
end

