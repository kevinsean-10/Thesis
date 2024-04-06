% classdef spiropt_fun
%     properties
%         m_cluster
%         k_cluster
%         m
%         k_max
%         gamma
%         epsilon
%         delta
%         radius
%         theta
%         dim
%         cluster_radius
%         cluster_center
%     end
% 
%     methods
%         function obj = spiropt_fun(m_cluster, k_cluster, m, k_max, gamma, epsilon, delta, radius, theta, dim)
%             obj.m_cluster = m_cluster;
%             obj.k_cluster = k_cluster;
%             obj.m = m;
%             obj.k_max = k_max;
%             obj.gamma = gamma;
%             obj.epsilon = epsilon;
%             obj.delta = delta;
%             obj.radius = radius;
%             obj.theta = theta;
%             obj.dim = dim;
%             obj.cluster_radius
%             obj.cluster_center
% 
%         end
%         function res = objective_function(obj,x)
%             f1 = exp(x(1)-x(2)) - sin(x(1)+x(2));
%             f2 = (x(1)*x(2))^2 - cos(x(1)+x(2));
%             F_array = [f1,f2];
%             res = 0;
%             for i = 1:length(F_array)
%                 res = res + abs(F_array(i));
%             end
%             res = -1 / (1 + res);
%         end
% 
% %         function f = objective_function(obj,x)
% %             f = 0;
% %             % Schwefel
% %             for i = 1:numel(x)
% %                 sum_sq = 0;
% %                 for j = 1:i
% %                     sum_sq = sum_sq + x(j);
% %                 end
% %                 f = f + sum_sq^2;
% %             end
% %             % 2^n Minima
% %             % for i = 1:numel(x)
% %             %     f = f + (x(i)^4 - 16*x(i)^2 + 5*x(i));
% %             % end
% %             % Rastrigin
% %             % for i = 1:numel(x)
% %             %     f = f + (x(i)^2 - 10*cos(2*pi*x(i)) + 10);
% %             % end
% %         end
% 
%         function z = generate_points(obj,nVar,N,LB,UB)
%             p = sobolset(nVar);%,'Skip',1e4,'Leap',1e2);
%             p = scramble(p,'MatousekAffineOwen');
%             A=net(p,N);
%             for i=1:nVar
%                 z(:,i)=round((LB(i)+(UB(i)-LB(i)).*A(:,i))*100)/100;
%             end
%         end
%         
%         function Rn_ij = generate_Rij(obj, i, j, dim, theta)
%             Rn_ij = eye(dim);
%             Rn_ij(i, i) = cos(theta);
%             Rn_ij(i, j) = -sin(theta);
%             Rn_ij(j, i) = sin(theta);
%             Rn_ij(j, j) = cos(theta);
%         end
% 
%         function Rn = generate_Rn(obj, dim, theta)
%             Rn = eye(dim);
%             for i = 1:dim-1
%                 for j=1:i
%                     product = eye(dim);
%                     product = product * obj.generate_Rij(dim-i, dim+1-j, dim, theta);
%                     Rn = Rn*product;
%                 end
%             end
%         end
% 
%         function new_set_of_points = update_point(obj,set_of_points, Sn)
%             fitness = zeros(size(set_of_points, 1), 1);
%             for i = 1:size(set_of_points, 1)
%                 fitness(i) = obj.objective_function(set_of_points(i, :));
%             end
%             [~, i_g] = min(fitness);
%             x_i_g = set_of_points(i_g, :);
%         
%             new_set_of_points = set_of_points;
%             dim = size(set_of_points, 2);
%         
%             for i = 1:size(new_set_of_points, 1)
%                 poin = Sn * set_of_points(i, :)' - (Sn - eye(dim)) * x_i_g';
%                 new_set_of_points(i, :) = poin';
%             end
%         end
% 
%         function err = iter_error(obj,old_set_of_points, new_set_of_points)
%             err = 0;
%             for i = 1:size(old_set_of_points, 1)
%                 diff = abs(norm(old_set_of_points(i, :)) - norm(new_set_of_points(i, :)));
%                 if diff > err
%                     err = diff;
%                 end
%             end
%         end
% 
%         function [return_points, return_points_value] = spiral_opt(obj, boundaries, npoint, theta, radius, max_iter, max_error, seed)
%             dim = size(boundaries, 1);
%             Rn = obj.generate_Rn(dim, theta);
%             Sn = radius * Rn;
%             iter = 0;
%             iter_points = cell(1, max_iter + 1);
%             iter_points{1} = obj.generate_points(dim, npoint, boundaries(:, 1), boundaries(:, 2));
%             while iter <= max_iter
%                 iter_points{iter + 2} = obj.update_point(iter_points{iter + 1}, Sn);
%                 error = obj.iter_error(iter_points{iter + 1}, iter_points{iter + 2});
%                 if error < max_error
%                     break;
%                 end
%                 iter = iter + 1;
%             end
%         
%             return_points = zeros(1, dim);
%             for i = 1:dim
%                 return_points(i) = mean(iter_points{iter + 1}(:, i));
%             end
%             return_points_value = obj.objective_function(return_points);
%         end
% 
% 
%         function [cluster_center, cluster_radius] = function_cluster(obj, y, cluster_center, cluster_radius)
%             dist_list = vecnorm(cluster_center - y, 2, 2);
%             [~, min_dist_id] = min(dist_list);
%             xc = cluster_center(min_dist_id, :);
%             xt = (xc + y) / 2;
%             
%             Fxt = obj.objective_function(xt);
%             Fxc = obj.objective_function(xc);
%             Fy = obj.objective_function(y);
%             
%             if (Fxt > Fy) && (Fxt > Fxc)
%                 cluster_center = [cluster_center; y];
%                 cluster_radius = [cluster_radius; vecnorm(y - xt)];
%             elseif (Fxt < Fy) && (Fxt < Fxc)
%                 cluster_center = [cluster_center; y];
%                 cluster_radius = [cluster_radius; vecnorm(y - xt)];
%                 [cluster_center, cluster_radius] = obj.function_cluster(xt, cluster_center, cluster_radius);
%             elseif Fy < Fxc
%                 cluster_center(min_dist_id, :) = y;
%             end
%             
%             cluster_radius(min_dist_id) = vecnorm(y - xt);
%         end
% 
%     end
% end
% 
% 

classdef Spiral_Optimization
    properties
        boundaries       
        m_cluster
        k_cluster
        m
        k_max
        radius
        theta
        gamma
        epsilon
        delta
        dim
        condition
        archive
        score
        final_root
        cluster_center
        cluster_radius
        iter_points
    end
    
    methods
        function obj = Spiral_Optimization(boundaries, m_cluster, k_cluster, m, k_max, radius, theta, gamma, epsilon, delta, seed)
            obj.boundaries = boundaries;
            obj.m_cluster = m_cluster;
            obj.k_cluster = k_cluster;
            obj.m = m;
            obj.k_max = k_max;
            obj.radius = radius;
            obj.theta = theta;
            obj.gamma = gamma;
            obj.epsilon = epsilon;
            obj.delta = delta;
            obj.dim = size(boundaries, 1);
            obj.condition = -1 + epsilon;
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

%         function points = generate_points(obj, npoint, low, high, sobol)
%             if nargin < 5
%                 sobol = true;
%             end
% 
%             if isequal(class(low), 'double')
%                 low = repmat(low, obj.dim, 1);
%                 high = repmat(high, obj.dim, 1);
%             end
% 
%             if sobol
%                 sampler = qmc.Sobol('d', obj.dim, 'scramble', true, 'seed', obj.seed);
%                 sample = sampler.random(npoint);
%                 points = qmc.scale('sample', sample, 'l_bounds', low, 'u_bounds', high);
%             else
%                 rng(obj.seed);
%                 points = rand(npoint, obj.dim);
%                 for i = 1:obj.dim
%                     points(:, i) = low(i) + (high(i) - low(i)) * points(:, i);
%                 end
%             end
%         end

        function points = generate_points(obj,npoint,LB,UB)
            p = sobolset(obj.dim);%,'Skip',1e4,'Leap',1e2);
            p = scramble(p,'MatousekAffineOwen');
            A = net(p,npoint);
            for i=1:obj.dim
                points(:,i)=round((LB(i)+(UB(i)-LB(i)).*A(:,i))*100)/100;
            end
        end
        
        function Rn_ij = generate_Rij(obj, i, j)
            Rn_ij = eye(obj.dim);
            Rn_ij(i, i) = cos(obj.theta);
            Rn_ij(i, j) = -sin(obj.theta);
            Rn_ij(j, i) = sin(obj.theta);
            Rn_ij(j, j) = cos(obj.theta);
        end
        
        function Rn = generate_Rn(obj)
            Rn = eye(obj.dim);
            for i = 1:obj.dim
                product = eye(obj.dim);
                for j = 1:i
                    product = product * obj.generate_Rij(obj.dim - i + 1, obj.dim + 1 - j);
                end
                Rn = Rn * product;
            end
        end

        function new_points = update_point(obj, set_of_points)
            Sn = obj.radius * obj.generate_Rn();
            fitness = zeros(size(set_of_points, 1), 1);
            for i = 1:size(set_of_points, 1)
                fitness(i) = obj.objective_function(set_of_points(i, :));
            end
            [~, i_g] = min(fitness);
            x_i_g = set_of_points(i_g, :);

            new_points = set_of_points;
            dim = size(set_of_points, 2);

            for i = 1:size(set_of_points, 1)
                poin = Sn * set_of_points(i, :)' - (Sn - eye(dim)) * x_i_g';
                new_points(i, :) = poin';
            end
        end

        function err = iter_error(obj, old_set_of_points, new_set_of_points)
            err = 0;
            for i = 1:size(old_set_of_points, 1)
                diff = abs(norm(old_set_of_points(i, :)) - norm(new_set_of_points(i, :)));
                if diff > err
                    err = diff;
                end
            end
        end

        function [return_points, return_points_value] = spiral_opt(obj, max_iter, max_error)
            iter = 0;
            cluster_points = cell(1, max_iter + 1);
            cluster_points{1} = obj.generate_points(obj.m_cluster, obj.boundaries(:, 1), obj.boundaries(:, 2), true);
            while iter <= max_iter
                cluster_points{iter + 2} = obj.update_point(cluster_points{iter + 1});
                error = obj.iter_error(cluster_points{iter + 1}, cluster_points{iter + 2});
                if error < max_error
                    break;
                end
                iter = iter + 1;
            end

            return_points = zeros(obj.dim, 1);
            for i = 1:obj.dim
                return_points(i) = mean(cluster_points{end}(:, i));
            end
            return_points_value = obj.objective_function(return_points');
        end

        function initialization(obj)
            obj.iter_points = cell(1, obj.k_cluster);
            obj.iter_points{1} = obj.generate_points(obj.m_cluster, obj.boundaries(:, 1), obj.boundaries(:, 2));
            fitness = zeros(obj.m_cluster, 1);
            for i = 1:obj.m_cluster
                fitness(i) = obj.objective_function(obj.iter_points{1}(i, :));
            end
            [~, best_idx] = min(fitness);
            x_prime = obj.iter_points{1}(best_idx, :);

            init_radius = (obj.boundaries(:, 2) - obj.boundaries(:, 1)) / 2;
            [~, id_rad] = min(init_radius);
            radii = init_radius(id_rad);

            obj.cluster_center = x_prime;
            obj.cluster_radius = radii;
        end
        function function_cluster(obj, y)
            dist_list = vecnorm(obj.cluster_center - y, 2, 2);
            [min_dist, min_dist_id] = min(dist_list);
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
        function clustering(obj)
            k = 0;
            while k < obj.k_cluster
                potential_cluster_center = [];
                F = zeros(obj.m_cluster, 1);
                for i = 1:obj.m_cluster
                    F(i) = obj.objective_function(obj.iter_points{k+1}(i, :)');
                end
                for i = 1:obj.m_cluster
                    exist_in_cluster_center = any(vecnorm(obj.iter_points{k+1}(i, :) - obj.cluster_center, 2, 2) < obj.epsilon);
                    if (F(i) < obj.gamma) && (~exist_in_cluster_center)
                        potential_cluster_center = [potential_cluster_center; obj.iter_points{k+1}(i, :)];
                    end
                end
                for i = 1:size(potential_cluster_center, 1)
                    obj.function_cluster(potential_cluster_center(i, :));
                end
                obj.iter_points{k+2} = obj.update_point(obj.iter_points{k+1});
                k = k + 1;
            end
        end
        function cluster_2Dvisualization(obj)
            if obj.dim ~= 2
                disp(['Dimension ', num2str(obj.dim), ' can be visualized using cluster_visualization2D.']);
                return;
            end

            figure;
            hold on;
            for i = 1:size(obj.cluster_center, 1)
                center = obj.cluster_center(i, :);
                radius = obj.cluster_radius(i);
                viscircles(center, radius, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1);
            end

            % Set axis limits
            xlim(obj.boundaries(1, :));
            ylim(obj.boundaries(2, :));

            % Add labels and title
            title('Cluster Visualization');
            xlabel('X-axis');
            ylabel('Y-axis');

            % Show the plot
            axis equal;
            grid on;
            hold off;
        end
        function final_root = root_elimination(obj, root_archive)
            if obj.dim == 1
                list_criteria = [root_archive{:}];
            else
                list_criteria = root_archive;
            end

            eligible_roots = [];
            for i = 1:numel(list_criteria)
                if obj.objective_function(list_criteria{i}) < obj.condition
                    eligible_roots = [eligible_roots; list_criteria{i}];
                end
            end

            id_duplicated_roots = [];
            for i = 1:numel(eligible_roots)
                for j = i+1:numel(eligible_roots)
                    if norm(eligible_roots{i} - eligible_roots{j}) < obj.delta
                        id_duplicated_roots = [id_duplicated_roots; i, j];
                    end
                end
            end
            id_duplicated_roots = unique(id_duplicated_roots, 'rows');

            deselected_id_duplicated_roots = [];
            for i = 1:size(id_duplicated_roots, 1)
                root_a = obj.objective_function(eligible_roots{id_duplicated_roots(i, 1)});
                root_b = obj.objective_function(eligible_roots{id_duplicated_roots(i, 2)});
                if root_a <= root_b
                    id_duplicated_root = id_duplicated_roots(i, 2);
                else
                    id_duplicated_root = id_duplicated_roots(i, 1);
                end
                deselected_id_duplicated_roots = [deselected_id_duplicated_roots; id_duplicated_root];
            end

            if ~isempty(deselected_id_duplicated_roots)
                unique_roots = true(1, numel(eligible_roots));
                unique_roots(deselected_id_duplicated_roots) = false;
                final_root = eligible_roots(unique_roots);
            else
                final_root = eligible_roots;
            end
        end

        function spiral_opt_evaluation(obj, verbose)
            if nargin < 2
                verbose = false;
            end
            
            obj.archive = {};
            obj.score = [];
            for i = 1:numel(obj.cluster_center)
                subbound = [obj.cluster_center{i}-obj.cluster_radius(i); obj.cluster_center{i}+obj.cluster_radius(i)];
                [root, root_score] = obj.spiral_opt(@obj.objective_function, subbound, obj.m, obj.k_max, obj.epsilon);
                obj.archive{i} = root;
                obj.score(i) = root_score;

                if verbose
                    fprintf('\n====== Cluster %d ======\n', i);
                    fprintf('Roots = \n');
                    disp(obj.archive);
                end
            end

            obj.final_root = obj.root_elimination(obj.archive);
        end
    end
end



