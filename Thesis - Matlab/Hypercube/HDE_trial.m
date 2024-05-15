clear;clc;

epsilon = 1e-5;
delta = 0.01;
pop_size = 100;
max_gen = 300;
dim = 2;
mutation_factor=0.1;
crossover_rate=0.5;
seed = 'shuffle';
print_stat = false;
verbose = true;
visual_properties = struct('show_visual',true, ...
    'save_visual', false, ...
    'file_name', 'hde.avi');

% how many parts/slices do you desire in each dimension?
parts = 200;

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);
% x1+10>=0; -x1+10>=0
% x2+10>=0; -x2+10>=0

hdeopt = HDE(boundaries,pop_size,parts,max_gen,mutation_factor,crossover_rate,epsilon,delta,seed);

[final_root,final_score] = hdeopt.DE_evaluation(verbose,print_stat)
hdeopt.visualization2D(visual_properties)
pause(5);
close;

%% Exporting Statistic
clear; clc;

epsilon = 1e-6;
delta = 0.01;
pop_size = 250;
parts = [100];
max_gen = [50,100,250];
dim = 2;
mutation_factor=0.1;
crossover_rate=0.5;
print_stat = false;
verbose = false;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'hde.avi');

boundaries = repmat([-10, 10], dim, 1);

max_iter = 100;

basename = 'hga';
us = '_';
extension = '.avi';
xlsx = '.xlsx';

disp('-start-')
sheet1 = []; 
sheet2 = [];
sheet3 = [];
for gm = 1: size(max_gen,2)
    for pt = 1: size(parts,2)
        disp("----------")
        fprintf('parts=%d\nmax_gen=%d\n',parts(pt),max_gen(gm))
        for iter=1:max_iter
            fprintf('Iteration: %d\n',iter)
            seed = 'shuffle';
            rngseed = rng(seed);
            print_seed = rngseed.Seed;
            visual_properties = struct('show_visual',false, ...
            'save_visual', false, ...
            'file_name', [basename,us,num2str(parts(pt)),us,num2str(max_gen(gm)),us,num2str(iter),us,num2str(rngseed.Seed),extension]);
            tic;

            hdeopt = HDE(boundaries,pop_size,parts(pt),max_gen(gm),mutation_factor,crossover_rate,epsilon,delta,seed);
            [final_root,final_score] = hdeopt.DE_evaluation(verbose,print_stat);

            if isempty(final_root)
                final_root = NaN(1,dim);
                final_score = NaN;
            end

            elapsed_time = toc;

            % 1st Sheet
            num_iter = iter*ones(size(final_root,1),1);
            sheet1 = [sheet1; num_iter,final_root,final_score'];

            % 2nd Sheet
            num_root = size(final_root,1);
            best_score = min(final_score);
            sheet2 = [sheet2; iter,num_root,best_score,elapsed_time];
            sheet3 = [sheet3; iter,print_seed];
            writematrix(sheet1 ,[basename,us,num2str(parts(pt)),us,num2str(max_gen(gm)),xlsx],'Sheet','Final Root')
            writematrix(sheet2,[basename,us,num2str(parts(pt)),us,num2str(max_gen(gm)),xlsx],'Sheet','Statistic')
            writematrix(sheet3,[basename,us,num2str(parts(pt)),us,num2str(max_gen(gm)),xlsx],'Sheet','Random Seed')
            close;
        end
    end
end

disp('-end-')