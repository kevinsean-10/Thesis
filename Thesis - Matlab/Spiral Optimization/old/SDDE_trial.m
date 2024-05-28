%% Trial
clear; clc;

dim = 2;
mutation_factor=0.8;
crossover_rate=0.1;
tau_d=0.4;
m_cluster = 250;
gamma = -0.2;
epsilon = 1e-6;
delta = 0.01;
k_cluster = 10;
m = 250;
k_max = 250;
seed = 'shuffle';
verbose = false;
print_stat = false;
visual_properties = struct('show_visual',false, ...
    'save_visual', false, ...
    'file_name', 'sdde.avi');

% Define boundaries
boundaries = repmat([-10, 10], dim, 1);

sddeopt = SDDE(boundaries,m_cluster,k_cluster,m,k_max, ...
                epsilon,delta,gamma,mutation_factor,crossover_rate,seed);

sddeopt.DE_evaluation(verbose,print_stat,visual_properties)

% %% Exporting Statistic
% clear; clc;
% 
% dim = 2;
% mutation_factor=0.8;
% crossover_rate=0.1;
% tau_d=0.4;
% m_cluster = [50,100,250];
% gamma = -0.2;
% epsilon = 1e-6;
% delta = 0.01;
% k_cluster = 10;
% m = 250;
% k_max = [50,100,250];
% boundaries = repmat([-10, 10], dim, 1);
% verbose = false;
% print_stat = false;
% 
% max_iter = 100;
% 
% basename = 'sdde';
% us = '_';
% extension = '.avi';
% xlsx = '.xlsx';
% 
% disp('-start-')
% sheet1 = []; sheet2 = []; sheet3 = [];
% for km = 1: size(k_max,2)
%     for mc = 1: size(m_cluster,2)
%         disp("----------")
%         fprintf('m_cluster=%d\nk_max=%d\n',m_cluster(mc),k_max(km))
%         for iter=1:max_iter
%             fprintf('Iteration: %d\n',iter)
%             seed = 'shuffle';
%             rngseed = rng(seed);
%             print_seed = rngseed.Seed;
%             visual_properties = struct('show_visual',false, ...
%             'save_visual', false, ...
%             'file_name', [basename,us,num2str(m_cluster(mc)),us,num2str(k_max(km)),us,num2str(iter),us,num2str(rngseed.Seed),extension]);
%             tic;
%             tic;
%             seed = 'shuffle';
% 
%             sddeopt = SDDE(boundaries,m_cluster(mc),k_cluster,m,k_max(km), ...
%                         epsilon,delta,gamma,mutation_factor,crossover_rate,seed);
%             [final_root,final_score] = sddeopt.DE_evaluation(verbose,print_stat,visual_properties);
% 
%             elapsed_time = toc;
%             if isempty(final_root)
%                 final_root = NaN(1,dim);
%                 final_score = NaN;
%             end
%             % 1st Sheet
%             num_iter = iter*ones(size(final_root,1),1);
%             sheet1 = [sheet1; num_iter,final_root,final_score'];
% 
%             % 2nd Sheet
%             num_root = size(final_root,1);
%             best_score = min(final_score);
%             sheet2 = [sheet2; iter,num_root,best_score,elapsed_time];
%             sheet3 = [sheet3; iter,print_seed];
%             writematrix(sheet1 ,[basename,us,num2str(m_cluster(mc)),us,num2str(k_max(km)),xlsx],'Sheet','Final Root')
%             writematrix(sheet2,[basename,us,num2str(m_cluster(mc)),us,num2str(k_max(km)),xlsx],'Sheet','Statistic')
%             writematrix(sheet3,[basename,us,num2str(m_cluster(mc)),us,num2str(k_max(km)),xlsx],'Sheet','Random Seed')
%             close;
%         end
%     end
% end
% 
% disp('-end-')
% 
% 
