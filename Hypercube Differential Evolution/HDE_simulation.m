%% Exporting Statistic
clear; clc;

squared_parts = [200,400];
parts = round(sqrt(squared_parts));
max_gen = [100,250,400];

epsilon = 1e-6;
tau_d = 0.1;
pop_size = 300;
dim = 3;
mutation_factor=0.831987414;
crossover_rate=0.887292888;
seed = 'shuffle';
print_stat = false;
verbose = false;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'hde.avi');

% Define boundaries
% % Problem 1
% boundaries = repmat([-10, 10], dim, 1);
% % Problem 2
% boundaries = [-1,3;-17,4];
% Problem 7
boundaries = [0,2;-10,10;-1,1];

max_iter = 110;

basename = 'hde';
us = '_';
extension = '.avi';
xlsx = '.xlsx';

disp('-start-')
failedIndices = [];
for gm = 1: size(max_gen,2)
    for pt = 1: size(parts,2)
        sheet1 = []; 
        sheet2 = [];
        sheet3 = [];
        disp("----------")
        fprintf('parts=%d\nsquared_parts=%d\nmax_gen=%d\n',parts(pt),squared_parts(pt),max_gen(gm))
        for iter=1:max_iter
            try
                fprintf('Iteration: %d\n',iter)
                seed = 'shuffle';
                rngseed = rng(seed);
                print_seed = rngseed.Seed;
                visual_properties = struct('show_visual',false, ...
                'save_visual', false, ...
                'file_name', [basename,us,num2str(squared_parts(pt)),us,num2str(max_gen(gm)),us,num2str(iter),us,num2str(rngseed.Seed),extension]);
                tic;

                hdeopt = HDE(boundaries,pop_size,parts(pt),max_gen(gm),mutation_factor,crossover_rate,epsilon,tau_d,seed);
                [final_root,final_score,num_cluster] = hdeopt.DE_evaluation(verbose,print_stat);
                num_root = size(final_root,1);
                best_score = min(final_score);

                if isempty(final_root)
                    final_root = NaN(1,dim);
                    final_score = NaN;
                end

                elapsed_time = toc;

                % 1st Sheet
                num_iter = iter*ones(size(final_root,1),1);
                sheet1 = [sheet1; num_iter,final_root,final_score'];

                % 2nd Sheet
                sheet2 = [sheet2; iter,num_root,best_score,num_cluster,elapsed_time];
                sheet3 = [sheet3; iter,print_seed];
                writematrix(sheet1 ,[basename,us,num2str(squared_parts(pt)),us,num2str(max_gen(gm)),xlsx],'Sheet','Final Root')
                writematrix(sheet2,[basename,us,num2str(squared_parts(pt)),us,num2str(max_gen(gm)),xlsx],'Sheet','Statistic')
                writematrix(sheet3,[basename,us,num2str(squared_parts(pt)),us,num2str(max_gen(gm)),xlsx],'Sheet','Random Seed')
                close;
            catch ME
                % If there is an error, display a warning and skip to the next iteration
                warning('Error untuk pt %d - gm %d - iter %d: %s', squared_parts(pt),max_gen(gm),iter, ME.message);
                % Store the index of the failed computation
                failedIndices = [failedIndices; [squared_parts(pt),max_gen(gm),iter]];
                continue;
            end
        end
    end
end
writematrix(failedIndices,['failed_indices',xlsx],'Sheet','Failed Indices')
disp('-end-')