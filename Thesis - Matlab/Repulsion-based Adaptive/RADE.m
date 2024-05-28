classdef RADE < handle
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
        function obj = RADE(boundaries,population_size, ...
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

        % % Problem 6
        % function F_array = system_equations(obj,x)
        %     f1 = x(1)^2 + x(1)*x(2) - 6;
        %     f2 = x(1)^2 + x(2)^3 + 2*x(1)*x(2)^2-3;
        %     F_array = [f1; f2];
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

        function [pop,BestSol] = generate_points(obj,npoint,boundaries,seed)
            rng(seed);
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
            F_array = repmat(F,VarSize);
            donor_vec = population(a).Position+F_array.*(population(b).Position-population(c).Position);
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

        function [subpop,closest_indices] = subpopulating(obj, individual, population, t)
            population_array = reshape([population.Position], obj.dim, [])';
            distances = sqrt(sum((population_array - individual) .^ 2, 2));
        
            % Get the indices of the individuals with the smallest distances
            [~, sorted_indices] = sort(distances);
            closest_indices = sorted_indices(1:t);
        
            % Form the subpopulation with the closest individuals
            subpop = population(closest_indices);
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

        function Fi = update_Fi(obj,hi)
            a = rand();
            Fi = 0.1*tan(pi * (a - 0.5)) + obj.memories_F(hi);
            if Fi > 1
                Fi = 1;
            elseif Fi <=0
                Fi = obj.update_Fi(hi);
            end
        end
        
        function [Fi, CRi] = update_parameter(obj)
            % OUTPUT
            % Fi: scaling factor
            % CRi: crossover rate

            % Randomly select an index
            hi = randi(obj.max_memories_size);
            % Generate Fi using the Cauchy distribution with the location parameter MF(hi) and scale 0.1
            Fi = obj.update_Fi(hi);
            % Generate CRi using the Gaussian distribution with mean MCR(hi) and standard deviation 0.1
            CRi = normrnd(obj.memories_CR(hi), 0.1);
            % Ensure CRi is within the range [0, 1]
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

        function [final_root,final_score] = DE_evaluation(obj,verbose,visual_properties)
            rng(obj.seed);
            [population,bestsol] = obj.generate_points(obj.population_size,obj.boundaries,obj.seed);
            
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

            k=1;
            updated_list = [0];
            memoriesF = [];
            memoriesCR = [];
            VarSize = [1 obj.dim];
            for gen = 1:obj.max_generation
                updated = 0;
                    
                for ind=1:obj.population_size
                    obj.archive = obj.update_archive(population(ind).Position,obj.archive);
                end
                
                S_F = [];
                S_CR = [];
                for i = 1:obj.population_size
                    [F_i, CR_i] = obj.update_parameter();
                    x_i = population(i).Position;
                    [subpop_i,~] = obj.subpopulating(x_i, population, obj.num_per_subpopulation);

                    permvec = randperm(obj.num_per_subpopulation); % Randomize the indices of the population
                    permvec(permvec == i) = [];
                    mutant = obj.mutate(subpop_i,F_i,VarSize,obj.boundaries, permvec);
                    trial = obj.crossover(x_i,mutant,CR_i);
                    trial_fitness = obj.fitness_function(trial);
                    [subpop_closest_trial,closest_indices] = obj.subpopulating(trial,population,2);
                    closest_trial = subpop_closest_trial(end);
                    closest_trial.Cost = obj.fitness_function(closest_trial.Position);
                    
                    if trial_fitness < closest_trial.Cost
                        population(closest_indices(end)).Position = trial;
                        population(closest_indices(end)).Cost = trial_fitness;
                        S_F = [S_F, F_i];
                        S_CR = [S_CR, CR_i];
                        updated = updated + 1;
                        if trial_fitness < bestsol.Cost
                            best_idx = i;
                            bestsol.Position = trial;
                            bestsol.Cost = trial_fitness;
                        end
                    end
                    if ~isempty(S_F) && ~isempty(S_CR)
                        obj.update_history(S_F, S_CR, k);
                        k = k + 1;
                        if k > obj.max_memories_size
                            k = 1;
                        end
                    end

                end

                if verbose
                    fprintf("=========Generation %d=========\n", gen);
                    disp("Archive:");
                    disp(obj.archive);
                end

                if visual_properties.show_visual == true || visual_properties.save_visual == true
                    xlim(obj.boundaries(1,:));
                    ylim(obj.boundaries(2,:));
                    % answ = [-6.437160, 0.155348;
                    %         -0.932122, 1.067870;
                    %         -0.155283, 6.439840;
                    %          0.163333, 6.122430;
                    %          0.667121, 0.690103;
                    %         -6.117110, -0.163476];
                    for i=1:size(XY_surf,1)
                        Z_surf(i) = obj.repulsion_function(XY_surf(i,:));
                    end
                    Z_surf = reshape(Z_surf, size(X_surf));
                    subplot(2, 2, 1);
                    surf(X_surf, Y_surf, Z_surf);
                    view(0, 0);
                    title("Tampak Samping")

                    subplot(2, 2, 2);
                    memoriesF = [memoriesF,obj.memories_F(k)];
                    memoriesCR = [memoriesCR,obj.memories_CR(k)];
                    hold on
                    plot(0:gen-1,memoriesF,'Color','magenta')
                    plot(0:gen-1,memoriesCR,'Color','green')
                    legend(["F",'CR'], 'Location', 'southwest')
                    hold off
                    ylim([0,1])
                    title("Parameter F dan CR")
                    grid on;

                    subplot(2, 2, 3);
                    plot(0:gen-1,updated_list,"Color",'red');
                    title("Tingkat Disrupsi")
                    grid on;
                    
                    subplot(2, 2, 4);
                    % scatter(answ(:,1),answ(:,2), 50,"MarkerEdgeColor","#EDB120","LineWidth",2);
                    % hold on;
                    population_array = reshape([population.Position], obj.dim, [])';
                    scatter(population_array(:,1), population_array(:,2), 5, 'filled','blue'); 
                    % hold off;
                    rectangle('Position',[obj.boundaries(:,1)',(obj.boundaries(:,2)-obj.boundaries(:,1))'],'EdgeColor','#FF0000')
                    grid on;
                end

                if visual_properties.show_visual == true || visual_properties.save_visual == true
                    % Adjust aspect ratio
                    % axis equal;
                    % pbaspect([diff(xlim()) diff(ylim()) 1]);
        
                    % Maximize figure window
                    set(gcf, 'WindowState', 'maximized');
                    if ~isempty(obj.archive)
                        hold on
                        plot(obj.archive(:,1), obj.archive(:,2),'*',"color",'#EDB120',"MarkerSize",50);
                        hold off
                    end
                    title(sprintf("Generation %d",gen))
                    pause(0.01)
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(writerObj, frame);
                    end
                end

                updated_list = [updated_list;updated];
            end
            final_root = obj.archive;
            final_score = zeros(1, size(final_root,1));
            for fin_iter = 1:size(final_root,1)
                final_score(fin_iter) = obj.objective_function(final_root(fin_iter, :));
            end
        end
    end
end