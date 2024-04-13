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
        max_archive_size
        memories_F
        memories_CR
        max_memories_size
        seed
    end
    
    methods
        function obj = RADE(boundaries,population_size,num_per_subpopulation,max_generation,max_archive_size,theta,tau_d,F_init,CR_init,max_memories_size,beta,rho,seed)
            obj.boundaries = boundaries;
            obj.dim = size(obj.boundaries,1);
            obj.population_size = population_size;
            obj.num_per_subpopulation = num_per_subpopulation;
            obj.max_archive_size = max_archive_size;
            obj.max_generation = max_generation;
            obj.theta = theta;
            obj.beta = beta;
            obj.rho = rho;
            obj.tau_d = tau_d;
            obj.max_memories_size = max_memories_size;
            obj.memories_F = ones(1, max_memories_size) * F_init;
            obj.memories_CR = ones(1, max_memories_size) * CR_init;
            obj.archive = [];
            obj.seed = seed;
        end

        function F_array = system_equations(obj,x)
            f1 = exp(x(1)-x(2)) - sin(x(1)+x(2));
            f2 = (x(1)*x(2))^2 - cos(x(1)+x(2));
            F_array = [f1; f2];
        end

%         function res = objective_function(obj, x)
%             F_array = obj.system_equations(x);
%             res = sum(abs(F_array));
%             res = -1 / (1 + res);
%         end

        function res = objective_function(obj,x)
            res = 0;
            F_array = obj.system_equations(x);
            for i = 1:numel(F_array)
                res = res + (F_array(i))^2;
            end
        end

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

        function dv_i = mutation_penalty(obj, x_i, subpop_i, boundaries, scaling_factor)
            % Mutation function with penalty for DE
            % Inputs:
            % x_i: target x_i
            % subpop_i: number of individuals closest to x_i
            % boundaries: boundaries/constraints of the function
            % scaling_factor: scaling factor of the function
            % Output:
            % dv_i: donor vector that has been mutated and penalized.
        
            % Generate three distinct individuals xr1, xr1, xr1 from the current population randomly
            pop_ids = 1:size(subpop_i, 1);
            indices_to_delete = find(all(subpop_i == x_i, 2)); % Ensure that x_i is excluded from the selected subpopulation
            subpop_ids_no_i = setdiff(pop_ids, indices_to_delete);
            subpop_i = subpop_i(subpop_ids_no_i, :);
        
            % Mutation form the donor/mutation vector
            dv_i = obj.mutate(subpop_i, scaling_factor);
        
            % Set penalty for every donor vector that violates the boundaries
            for j = 1:size(dv_i, 2)
                if dv_i(j) < boundaries(j, 1)
                    dv_i(j) = (x_i(j) + boundaries(j, 1)) / 2;
                elseif dv_i(j) > boundaries(j, 2)
                    dv_i(j) = (x_i(j) + boundaries(j, 2)) / 2;
                end
            end
        end

        function subpop = subpopulating(obj, individual, population, t)
            % Calculate Euclidean distances and select t closest individuals.
            % Inputs:
            % individual: individual target
            % population: population where the target is
            % t: max number of units in a subpopulation
            % return_index: optional, flag to return indices of closest individuals
            % show_distances: optional, flag to print distances
            % Outputs:
            % closest_indices: array of closest index of the target individual
            % subpop: array of closest individuals of the target individual
        
            % Calculate the Euclidean distances from the individual to all others in the population
            distances = sqrt(sum((population - individual) .^ 2, 2));
        
            % Get the indices of the individuals with the smallest distances
            [~, sorted_indices] = sort(distances);
            closest_indices = sorted_indices(1:t);
        
            % Form the subpopulation with the closest individuals
            subpop = population(closest_indices, :);
        end

        function archive = update_archive(obj, x)
            % Update the archive with a new individual x
            % Inputs:
            % x: individual
            % Outputs:
            % archive: updated archive
        
            f_x = obj.objective_function(x)
            s = size(obj.archive,1); % current archive size

            if f_x < obj.theta % x is a root
                if s == 0 % archive is empty
                    obj.archive = [obj.archive; x];
                    s = s + 1;
                else
                    % Find the closest solution x_prime in the archive to x in the decision space
                    dist_list = vecnorm(x - obj.archive, 2, 2);
                    [dist_min, idx_min] = min(dist_list);
                    x_prime = obj.archive(idx_min, :);
                    f_x_prime = obj.objective_function(x_prime);
        
                    if dist_min < obj.tau_d % x and x_prime are too close
                        if f_x < f_x_prime
                            x_prime = x;
                            obj.archive(idx_min, :) = x_prime;
                        end
                    else
                        if s < obj.max_archive_size
                            obj.archive = [obj.archive; x];
                            s = s + 1;
                        else % archive is full
                            if f_x < f_x_prime
                                x_prime = x;
                                obj.archive(idx_min, :) = x_prime;
                            end
                        end
                    end
                end
            end
        
            archive = obj.archive;
        end
        
        function [Fi, CRi] = update_parameter(obj)
            % OUTPUT
            % Fi: scaling factor
            % CRi: crossover rate

            % Randomly select an index
            hi = randi(obj.max_memories_size);
            % Generate Fi using the Cauchy distribution with the location parameter MF(hi) and scale 0.1
%             Fi = tan(pi * (rand - 0.5)) + obj.memories_F(hi)
            Fi = obj.memories_F(hi) + 0.1*trnd(1,1);
            % Generate CRi using the Gaussian distribution with mean MCR(hi) and standard deviation 0.1
            CRi = normrnd(obj.memories_CR(hi), 0.1);
            % Ensure CRi is within the range [0, 1] and Fi is within the range [0,2]
            Fi = min(max(Fi, 0), 1);
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
            end
            result = numerator / denominator;

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


    end
end

