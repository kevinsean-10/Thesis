%% Trial
clear; clc; close all

% Problem 1 Parameter
pop_size=250;
max_gen=250;
F_init=0.5;
CR_init=0.5;
seed = 'shuffle';
num_l=10;
theta=1e-6;
tau_d=0.4;
s_max=20;
print_gen=false;
Hm = 50;
dim = 2;
beta = 1;
rho = 0.1;
boundaries = repmat([-10, 10], dim, 1);

% % Problem 2 Parameter
% pop_size=2000;
% max_gen=300;
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
% boundaries = [-1,3;-17,4];

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

rngseed = rng(seed);
disp(rngseed.Seed)
basename = 'rade';
us = '_';
extension = '.avi';

visual_properties = struct('show_visual',true, ...
                            'save_visual', false, ...
                            'file_name', [basename,us,num2str(pop_size),us,num2str(max_gen),us,num2str(rngseed.Seed),extension]);

% Problem 1 & 2
radeopt = RADE(boundaries,pop_size,num_l,max_gen,s_max,theta,tau_d,F_init,CR_init,Hm,beta,rho,seed);

% % Problem 3
% radeopt = RADE3(boundaries,pop_size,num_l,max_gen,s_max,theta,tau_d,F_init,CR_init,Hm,beta,rho,seed);

[final_root,final_score] = radeopt.DE_evaluation(print_gen,visual_properties)

% pause(2);
% close;

