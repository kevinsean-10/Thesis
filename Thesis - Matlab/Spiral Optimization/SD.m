classdef SD < handle
    properties
        boundaries
        dim
        seed
        epsilon
        gamma
        theta
        radius
        delta
        m_cluster
        k_cluster
        m
        k_max
        cluster_radius
        cluster_center
        cluster_iter_points
    end
    
    methods
        function obj = SD(boundaries,m_cluster,k_cluster,m,k_max,epsilon,delta,gamma,theta,radius,seed)
            obj.boundaries = boundaries;
            obj.dim = size(obj.boundaries,1);
            obj.m_cluster = m_cluster;
            obj.k_cluster = k_cluster;
            obj.m = m;
            obj.k_max = k_max;
            obj.epsilon = epsilon;
            obj.delta = delta;
            obj.gamma = gamma;
            obj.theta = theta;
            obj.seed = seed;
            obj.radius = radius;
            obj = initialization(obj);
        end
        
        function F_array = system_equations(obj,x)
            f1 = exp(x(1)-x(2)) - sin(x(1)+x(2));
            f2 = x(1)^2*x(2)^2 - cos(x(1)+x(2));
            F_array = [f1; f2];
        end

        function res = objective_function(obj, x)
            F_array = obj.system_equations(x);
            res = sum(abs(F_array));
            res = -1 / (1 + res);
        end

        % FUNGSI BARRIER
        function Phi = phi_func(obj, point)
            h_lowerbound = @(x) x' - obj.boundaries(:,1);
            h_upperbound = @(x) -(x' - obj.boundaries(:,2));
            h_lower = h_lowerbound(point);
            h_upper = h_upperbound(point);
            Phi = 0;
            for i=1:obj.dim
                Phi = Phi + log(h_lower(i)) + log(h_upper(i));
            end
            Phi = -Phi;
        end

        function Beta = beta_func(obj,point,mu)
            Beta = obj.objective_function(point) + mu*obj.phi_func(point);
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

        function eligible_points = inside_boundaries(obj,points, boundaries)
            criterion = all(points >= boundaries(:, 1)' & points <= boundaries(:, 2)',2);
            eligible_points = points(criterion,:);
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
            for i = 1:obj.dim-1
                product = eye(obj.dim);
                for j = 1:i
                    product = product * obj.generate_Rij(obj.dim-i, obj.dim+1-j);
                end
                Rn = Rn * product;
            end
        end

        function [new_set_of_points,x_i_g] = update_point(obj, set_of_points)
            Sn = obj.radius * obj.generate_Rn();
            fitness = zeros(size(set_of_points, 1), 1);
            for i = 1:size(set_of_points, 1)
                fitness(i) = obj.objective_function(set_of_points(i, :));
            end
            [~, i_g] = min(fitness);
            x_i_g = set_of_points(i_g, :);
        
            new_set_of_points = set_of_points;
            dimension = size(set_of_points, 2);
        
            for i = 1:size(new_set_of_points, 1)
                poin = Sn * set_of_points(i, :)' - (Sn - eye(dimension)) * x_i_g';
                new_set_of_points(i, :) = poin';
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

        function [return_points, return_points_value] = spiral_opt(obj, boundaries, n_point, max_iter, max_error,seed)
            points = obj.generate_points(n_point,boundaries,seed);
            iter = 0;
            while iter <= max_iter 
                [new_points,x_star] = obj.update_point(points);
                error = obj.iter_error(points, new_points);
                if error < max_error
                    break;
                end
                iter = iter + 1;
                points = new_points;
            end
            return_points = x_star;
            return_points_value = obj.objective_function(return_points);
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

        function clustering(obj, visual_properties)
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
                    exist_in_cluster_center = any(vecnorm(obj.cluster_iter_points(i,:) - obj.cluster_center) < obj.delta);
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
                        plot(obj.cluster_center(c,1), obj.cluster_center(c,2), '*',Color='magenta');
                        viscircles(obj.cluster_center(c,:), obj.cluster_radius(c), Color="#00FFFF", Linestyle='-.');
                        hold off
                    end
                    pause(0.25)
                    
                    if visual_properties.save_visual == true
                        frame = getframe(gcf);
                        writeVideo(writerObj, frame);
                    end
                end
        
                [obj.cluster_iter_points] = obj.update_point(obj.cluster_iter_points);
        
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

        function [final_root,final_score] = spiral_opt_evaluation(obj,verbose,visual_properties)
            obj.initialization();
            obj.clustering(visual_properties);
            if verbose == true
                fprintf("%d roots found!\n",numel(obj.cluster_radius))
            end
            
            archive = [];
            score = [];
            
            for i = 1:length(obj.cluster_center)
                subbound = [obj.cluster_center(i,:)' - obj.cluster_radius(i), obj.cluster_center(i,:)' + obj.cluster_radius(i)];
                [root, root_score] = obj.spiral_opt(subbound, obj.m, obj.k_max, obj.epsilon,obj.seed);
                inside_root = obj.inside_boundaries(root,obj.boundaries);
                if ~isempty(inside_root)
                    archive = [archive;root];
                    score = [score;root_score];
                end
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