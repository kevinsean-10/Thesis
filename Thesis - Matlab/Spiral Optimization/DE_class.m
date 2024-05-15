classdef DE_class < handle
    properties
        dim
        boundaries
        max_iter
        population_size
        scaling_factor
        crossover_prob
        seed
    end
    
    methods
        function obj = DE_class(boundaries,max_iter,population_size, ...
                scaling_factor,crossover_prob,seed)
            obj.boundaries = boundaries;
            obj.dim = size(boundaries,1);
            obj.max_iter = max_iter;
            obj.population_size = population_size;
            obj.scaling_factor = scaling_factor;
            obj.crossover_prob = crossover_prob;
            obj.seed = seed;
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

    end
end

