clearvars;
close all;

save = false;
draw=false; %Does not work with uniform & genie files
verbose=true;

%------------------------Input parameters------------------------%
parameters.J=10; %number of jobs
parameters.lambda_J=0.001; %arrival rate of jobs, in [job/slots]
parameters.feedback_interval=round(1./parameters.lambda_J); %time interval between two feedbacks, in [slots]

parameters.Omega_forward = 1.5; %a priori redundancy ratio
parameters.Omega = parameters.Omega_forward; %redundancy ratio
parameters.N = 100; %size of the matrix
parameters.m=50; %number of partitions
parameters.s=40;
parameters.t=parameters.m/parameters.s;
parameters.PP=1000; %number of available workers
parameters.Ns=100; %number of possible (s,t)
parameters.mu_rate_vec_init=1000*rand(1,parameters.PP); %service rate of the workers, in [comp/slots]
parameters.var_vec_init=2./parameters.mu_rate_vec_init.^2;
parameters.mu_enc = 10000; %encoding service rate 
parameters.mu_dec = 10000; %decoding service rate
parameters.c_rate_vec_init = 200*rand(1,parameters.PP); %communication rate of the workers, in [bits/slots]
parameters.Theta =2; %margin for workers' initial choice
parameters.beta = 0.01; %for update of mu
parameters.kappa = 0.2; %lower bound on utilization, sum phi_l = kappa
%parameters.fact_time_out = 2;
%% simulation
simu_results = main_solution_queue_fct_with_purging(parameters,save,draw,verbose);
