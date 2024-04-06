classdef func_cluster    
    properties
        cluster_center
        cluster_radius
    end
    
    methods
        function obj = func_cluster(cluster_center,cluster_radius)
            obj.cluster_center = cluster_center;
            obj.cluster_radius = cluster_radius;
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
    end
end

