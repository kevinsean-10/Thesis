%% Exporting Statistic
clear; clc;

% % Problem 1 Parameter
% pop_size=250;
% max_gen=250;
% F_init=0.5;
% CR_init=0.5;
% seed = 'shuffle';
% num_l=10;
% theta=1e-6;
% tau_d=0.4;
% s_max=20;
% print_gen=false;
% Hm = 50;
% dim = 2;
% beta = 1;
% rho = 0.1;
% boundaries = repmat([-10, 10], dim, 1);

% Problem 2 Parameter
pop_size=2000;
max_gen=300;
dim = 2;
mutation_factor=0.871070813;
crossover_rate=0.853720936;
tau_d=0.1;
gamma = -0.3;
epsilon = 1e-6;
delta = 0.01;
k_cluster = 10;
m = 300;
verbose = false;
print_stat = false;
seed = 'shuffle';
boundaries = [-1,3;-17,4];

% % Problem 3 Parameter
% pop_size=400;
% max_gen=400;
% F_init=0.5;
% CR_init=0.5;
% num_l=10;
% theta=1e-6;
% tau_d=0.1;
% s_max=20;
% print_gen=false;
% Hm = 50;
% dim = 2;
% beta = 1;
% seed = 'shuffle';
% rho = 0.1;
% boundaries = [0,2;-10,10;-1,1];

max_iter = 20;

basename = 'rade';
us = '_';
extension = '.avi';
xlsx = '.xlsx';

disp('-start-')
for ps = 1: size(pop_size,2)
    for mg = 1: size(max_gen,2)
        sheet1 = []; 
        sheet2 = [];
        sheet3 = [];
        disp("----------")
        fprintf('pop_size=%d\nmax_gen=%d\n',pop_size(ps),max_gen(mg))
        for iter=1:max_iter
            fprintf('Iteration: %d\n',iter)
            seed = 'shuffle';
            rngseed = rng(seed);
            print_seed = rngseed.Seed;
            visual_properties = struct('show_visual',false, ...
            'save_visual', false, ...
            'file_name', [basename,us,num2str(pop_size(ps)),us,num2str(max_gen(mg)),us,num2str(iter),us,num2str(rngseed.Seed),extension]);
            tic;

            radeopt = RADE(boundaries,pop_size(ps),num_l,max_gen(mg),s_max,theta,tau_d,F_init,CR_init,Hm,beta,rho,seed);
            [final_root,final_score] = radeopt.DE_evaluation(print_gen,visual_properties);

            elapsed_time = toc;

            if isempty(final_root)
                final_root = NaN(1,dim);
                final_score = NaN;
            end

            mean_F = mean(radeopt.memories_F);
            std_F = std(radeopt.memories_F);
            mean_CR = mean(radeopt.memories_CR);
            std_CR = std(radeopt.memories_CR);

            % 1st Sheet
            num_iter = iter*ones(size(final_root,1),1);
            sheet1 = [sheet1; num_iter,final_root,final_score'];

            % 2nd Sheet
            num_root = size(final_root,1);
            best_score = min(final_score);
            sheet2 = [sheet2; iter,num_root,best_score,mean_F,std_F,mean_CR,std_CR,elapsed_time];
            sheet3 = [sheet3; iter,print_seed];
            writematrix(sheet1 ,[basename,us,num2str(pop_size(ps)),us,num2str(max_gen(mg)),xlsx],'Sheet','Final Root')
            writematrix(sheet2,[basename,us,num2str(pop_size(ps)),us,num2str(max_gen(mg)),xlsx],'Sheet','Statistic')
            writematrix(sheet3,[basename,us,num2str(pop_size(ps)),us,num2str(max_gen(mg)),xlsx],'Sheet','Random Seed')
            close;
        end
    end
end

disp('-end-')