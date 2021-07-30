clearvars;
close all;

saved = false;
draw=false; %should not be set to true, otherwise too many figures.
verbose=false;
%------------------------Input parameters------------------------%
parameters.J=10; %number of jobs
parameters.lambda_J=0.001; %arrival rate of jobs, in [job/slots]
parameters.N = 100; %size of the matrix
parameters.m=50; %number of partitions
parameters.s=40; 
parameters.t=parameters.m/parameters.s;
parameters.PP=1000; %number of available workers
parameters.Ns=100; %number of possible (s,t)
parameters.mu_enc = 100000; %encoding service rate 
parameters.mu_dec = 1000000; %decoding service rate
parameters.Theta =2; %margin for workers' initial choice

parameters.beta = 0.01; %for update of mu
parameters.kappa = 0.2; %lower bound on utilization, sum phi_l = kappa
%parameters.fact_time_out = 2; %number of standard deviation for time_out
%parameters.Omega_forward = 1.05; %a priori redundancy ratio

NOmega = 3;
Omega_vec = linspace(1,1.7,NOmega);
Niter = 5; %for average

Dmean_vec_uni_no_purging = zeros(Niter,NOmega);
Dmax_vec_uni_no_purging = zeros(Niter,NOmega);
Thr_vec_uni_no_purging = zeros(Niter,NOmega);
glo_comm_load_ratio_vec_uni_no_purging = zeros(Niter,NOmega);
glo_comp_load_ratio_vec_uni_no_purging = zeros(Niter,NOmega);
idle_vec_uni_no_purging = zeros(Niter,NOmega);

Dmean_vec_uni_with_purging = zeros(Niter,NOmega);
Dmax_vec_uni_with_purging = zeros(Niter,NOmega);
Thr_vec_uni_with_purging = zeros(Niter,NOmega);
glo_comm_load_ratio_vec_uni_with_purging = zeros(Niter,NOmega);
glo_comp_load_ratio_vec_uni_with_purging = zeros(Niter,NOmega);
idle_vec_uni_with_purging = zeros(Niter,NOmega);

Dmean_vec = zeros(Niter,NOmega);
Dmax_vec = zeros(Niter,NOmega);
Thr_vec = zeros(Niter,NOmega);
glo_comm_load_ratio_vec = zeros(Niter,NOmega);
glo_comp_load_ratio_vec = zeros(Niter,NOmega);
idle_vec = zeros(Niter,NOmega);

Dmean_vec_no_purging = zeros(Niter,NOmega);
Dmax_vec_no_purging = zeros(Niter,NOmega);
Thr_vec_no_purging = zeros(Niter,NOmega);
glo_comm_load_ratio_vec_no_purging = zeros(Niter,NOmega);
glo_comp_load_ratio_vec_no_purging = zeros(Niter,NOmega);
idle_vec_no_purging = zeros(Niter,NOmega);

Dmean_vec_gini = zeros(Niter,NOmega);
Dmax_vec_gini = zeros(Niter,NOmega);
Thr_vec_gini = zeros(Niter,NOmega);
glo_comm_load_ratio_vec_gini = zeros(Niter,NOmega);
glo_comp_load_ratio_vec_gini = zeros(Niter,NOmega);
idle_vec_gini = zeros(Niter,NOmega);

Dmean_vec_gini_no_purging = zeros(Niter,NOmega);
Dmax_vec_gini_no_purging = zeros(Niter,NOmega);
Thr_vec_gini_no_purging = zeros(Niter,NOmega);
glo_comm_load_ratio_vec_gini_no_purging = zeros(Niter,NOmega);
glo_comp_load_ratio_vec_gini_no_purging = zeros(Niter,NOmega);
idle_vec_gini_no_purging = zeros(Niter,NOmega);


time_vec = zeros(Niter,NOmega);


    
fig1 = figure(1);  
fig2 = figure(2);  
fig3 = figure(3); 
fig4 = figure(4);
fig5 = figure(5); 
fig6 = figure(6); 
fig7 = figure(7); 
fig8 = figure(8); 
for niter = 1:Niter
    parameters.mu_rate_vec_init=2500*rand(1,parameters.PP); %service rate of the workers, in [comp/slots]
    parameters.c_rate_vec_init = 1000*rand(1,parameters.PP); %communication rate of the workers, in [bits/slots]
    for nn = 1:NOmega
        parameters.Omega_forward = Omega_vec(nn); %redundancy factor
        parameters.Omega = parameters.Omega_forward;
        
        tic;
        simu_results = main_uniform_queue_fct_no_purging(parameters,saved,draw,verbose);
        Dmean_vec_uni_no_purging(niter,nn)=simu_results.analysis.mean_total_time;
        Dmax_vec_uni_no_purging(niter,nn)=simu_results.analysis.max_total_time;
        Thr_vec_uni_no_purging(niter,nn) = Omega_vec(nn);
        glo_comp_load_ratio_vec_uni_no_purging(niter,nn) = simu_results.analysis.glo_comp_load_ratio;
        %idle_vec(niter,nn) = simu_results.analysis.percentage_idle_tot;
        
        simu_results = main_uniform_queue_fct_with_purging(parameters,saved,draw,verbose);
        Dmean_vec_uni_with_purging(niter,nn)=simu_results.analysis.mean_total_time;
        Dmax_vec_uni_with_purging(niter,nn)=simu_results.analysis.max_total_time;
        Thr_vec_uni_with_purging(niter,nn) = Omega_vec(nn);
        glo_comp_load_ratio_vec_uni_with_purging(niter,nn) = simu_results.analysis.glo_comp_load_ratio;
        %idle_vec(niter,nn) = simu_results.analysis.percentage_idle_tot;
        
        simu_results = main_solution_queue_fct_with_purging(parameters,saved,draw,verbose);
        Dmean_vec(niter,nn)=simu_results.analysis.mean_total_time;
        Dmax_vec(niter,nn)=simu_results.analysis.max_total_time;
        Thr_vec(niter,nn) = Omega_vec(nn);
        glo_comp_load_ratio_vec(niter,nn) = simu_results.analysis.glo_comp_load_ratio;
        %idle_vec(niter,nn) = simu_results.analysis.percentage_idle_tot;
        
        simu_results = main_solution_queue_fct_no_purging(parameters,saved,draw,verbose);
        Dmean_vec_no_purging(niter,nn)=simu_results.analysis.mean_total_time;
        Dmax_vec_no_purging(niter,nn)=simu_results.analysis.max_total_time;
        Thr_vec_no_purging(niter,nn) = Omega_vec(nn);
        glo_comp_load_ratio_vec_no_purging(niter,nn) = simu_results.analysis.glo_comp_load_ratio;
        %idle_vec(niter,nn) = simu_results.analysis.percentage_idle_tot;
        
        simu_results = main_solution_queue_fct_gini_with_purging(parameters,saved,draw,verbose);
        Dmean_vec_gini(niter,nn)=simu_results.analysis.mean_total_time;
        Dmax_vec_gini(niter,nn)=simu_results.analysis.max_total_time;
        Thr_vec_gini(niter,nn) = Omega_vec(nn);
        glo_comp_load_ratio_vec_gini(niter,nn) = simu_results.analysis.glo_comp_load_ratio;
        %idle_vec(niter,nn) = simu_results.analysis.percentage_idle_tot;
        
        simu_results = main_solution_queue_fct_gini_no_purging(parameters,saved,draw,verbose);
        Dmean_vec_gini_no_purging(niter,nn)=simu_results.analysis.mean_total_time;
        Dmax_vec_gini_no_purging(niter,nn)=simu_results.analysis.max_total_time;
        Thr_vec_gini_no_purging(niter,nn) = Omega_vec(nn);
        glo_comp_load_ratio_vec_gini_no_purging(niter,nn) = simu_results.analysis.glo_comp_load_ratio;
        %idle_vec(niter,nn) = simu_results.analysis.percentage_idle_tot;
        
        time = toc;
        time_vec(niter,nn)=time;
        
        mean_time = sum(sum(time_vec))/((niter-1)*NOmega+nn);
        remaining_step = (Niter-niter)*NOmega+NOmega-nn;
        fprintf('Iter %d/%d, Omega %d/%d, time %4.0f s \n',niter, Niter,nn, NOmega,time);
        fprintf('Remaining time: +/- %4.0f s \n', remaining_step*mean_time);
    end
    fprintf('========================================\n');
    close(fig3)
    fig3 = figure(3); 
    mean_delay_uni_no_purging = mean(Dmean_vec_uni_no_purging(1:niter,:),1);
    mean_thr_uni_no_purging = mean(glo_comp_load_ratio_vec_uni_no_purging(1:niter,:),1);
    scatter(reshape(Dmean_vec_uni_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
    hold on;
    plot(mean_delay_uni_no_purging,mean_thr_uni_no_purging,'o-r','markerfacecolor','r');
    %scatter(mean_delay_uni_no_purging(1),mean_thr_uni_no_purging(1),500,'r')
    title('Computational load vs Delay - uniform without purging');
    xlabel('Mean delay');
    ylabel('Computational load (number of tasks/optimum number)');
    
    close(fig4)
    fig4 = figure(4); 
    mean_delay_uni_with_purging = mean(Dmean_vec_uni_with_purging(1:niter,:),1);
    mean_thr_uni_with_purging = mean(glo_comp_load_ratio_vec_uni_with_purging(1:niter,:),1);
    scatter(reshape(Dmean_vec_uni_with_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_with_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
    hold on;
    plot(mean_delay_uni_with_purging,mean_thr_uni_with_purging,'o-k','markerfacecolor','k');
    %scatter(mean_delay_uni_with_purging(1),mean_thr_uni_with_purging(1),500,'r')
    title('Computational load vs Delay - uniform with purging');
    xlabel('Mean delay');
    ylabel('Computational load (number of tasks/optimum number)');
    
    
    close(fig5)
    fig5 = figure(5);
    mean_delay_no_purging = mean(Dmean_vec_no_purging(1:niter,:),1);
    mean_thr_no_purging = mean(glo_comp_load_ratio_vec_no_purging(1:niter,:),1);
    scatter(reshape(Dmean_vec_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
    hold on;
    plot(mean_delay_no_purging,mean_thr_no_purging,'o-c','markerfacecolor','c');
   % scatter(mean_delay(1),mean_thr(1),500,'r')
    title('Computational load vs Delay - optimum without purging');
    xlabel('Mean delay');
    ylabel('Computational load (number of tasks/optimum number)');
    
    
    close(fig6)
    fig6 = figure(6); 
    mean_delay = mean(Dmean_vec(1:niter,:),1);
    mean_thr = mean(glo_comp_load_ratio_vec(1:niter,:),1);
    scatter(reshape(Dmean_vec(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
    hold on;
    plot(mean_delay,mean_thr,'o-b','markerfacecolor','b');
   % scatter(mean_delay(1),mean_thr(1),500,'r')
    title('Computational load vs Delay - optimum with purging');
    xlabel('Mean delay');
    ylabel('Computational load (number of tasks/optimum number)');
    
    close(fig2)
    fig2 = figure(2); 
    mean_delay_gini = mean(Dmean_vec_gini(1:niter,:),1);
    mean_thr_gini = mean(glo_comp_load_ratio_vec_gini(1:niter,:),1);
    scatter(reshape(Dmean_vec_gini(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
    hold on;
    plot(mean_delay_gini,mean_thr_gini,'o-m','markerfacecolor','m');
   % scatter(mean_delay(1),mean_thr(1),500,'r')
    title('Computational load vs Delay - gini with purging');
    xlabel('Mean delay');
    ylabel('Computational load (number of tasks/optimum number)');
    
    close(fig1)
    fig1 = figure(1); 
    mean_delay_gini_no_purging = mean(Dmean_vec_gini_no_purging(1:niter,:),1);
    mean_thr_gini_no_purging = mean(glo_comp_load_ratio_vec_gini_no_purging(1:niter,:),1);
    scatter(reshape(Dmean_vec_gini_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
    hold on;
    plot(mean_delay_gini_no_purging,mean_thr_gini_no_purging,'o-g','markerfacecolor','g');
   % scatter(mean_delay(1),mean_thr(1),500,'r')
    title('Computational load vs Delay - gini without purging');
    xlabel('Mean delay');
    ylabel('Computational load (number of tasks/optimum number)');
    
    close(fig7)
    fig7 = figure(7); 
    plot(mean_delay_gini_no_purging,mean_thr_gini_no_purging,'o-g','markerfacecolor','g');
    hold on;
    plot(mean_delay_gini,mean_thr_gini,'o-m','markerfacecolor','m');
    plot(mean_delay_uni_no_purging,mean_thr_uni_no_purging,'o-r','markerfacecolor','r');
    plot(mean_delay_uni_with_purging,mean_thr_uni_with_purging,'o-k','markerfacecolor','k');
    plot(mean_delay_no_purging,mean_thr_no_purging,'o-c','markerfacecolor','c');
    plot(mean_delay,mean_thr,'o-b','markerfacecolor','b');
%scatter(mean_delay(1),mean_thr(1),500,'r')
    title('Computational load vs Delay');
    xlabel('Mean delay');
    ylabel('Computational load (number of tasks/optimum number)');
    legend('gini without purging', 'gini with purging','uni without purging', 'uni with purging', 'our without purging', 'our with purging');
    
    close(fig8)
    fig8 = figure(8); 
    scatter(reshape(Dmean_vec_gini_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
    hold on;
    scatter(reshape(Dmean_vec_gini(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini(1:niter,:),1,niter*NOmega),10,'m'); %to see the fog
    scatter(reshape(Dmean_vec_uni_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_no_purging(1:niter,:),1,niter*NOmega),10,'r'); %to see the fog
    scatter(reshape(Dmean_vec_uni_with_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_with_purging(1:niter,:),1,niter*NOmega),10,'k'); %to see the fog
    scatter(reshape(Dmean_vec_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_no_purging(1:niter,:),1,niter*NOmega),10,'c'); %to see the fog
    scatter(reshape(Dmean_vec(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec(1:niter,:),1,niter*NOmega),10,'b'); %to see the fog
    %scatter(mean_delay(1),mean_thr(1),500,'r')
    title('Computational load vs Delay');
    xlabel('Mean delay');
    ylabel('Computational load (number of tasks/optimum number)');
    legend('gini without purging', 'gini with purging', 'uni without purging', 'uni with purging', 'our without purging', 'our with purging');

    drawnow
end

%% save information
path = pwd;
path = strcat(path,'\simu_results');
c = num2cell(clock);
name = ['Distributed_computing_trdoff_simu2',...
                num2str(c{1}),'_', num2str(c{2}),'_', num2str(c{3}),'_', num2str(c{4}),'_', num2str(c{5})];
save([path, '/', name, '.mat']);
% 

 %% 
close(fig3)
fig3 = figure(3); 
scatter(reshape(Dmean_vec_uni_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
hold on;
plot(mean_delay_uni_no_purging,mean_thr_uni_no_purging,'o-r','markerfacecolor','r');
%scatter(mean_delay_uni_no_purging(1),mean_thr_uni_no_purging(1),500,'r')
title('Computational load vs Delay - uniform without purging');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');

close(fig4)
fig4 = figure(4); 
scatter(reshape(Dmean_vec_uni_with_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_with_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
hold on;
plot(mean_delay_uni_with_purging,mean_thr_uni_with_purging,'o-k','markerfacecolor','k');
%scatter(mean_delay_uni_with_purging(1),mean_thr_uni_with_purging(1),500,'r')
title('Computational load vs Delay - uniform with purging');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');


close(fig5)
fig5 = figure(5);
scatter(reshape(Dmean_vec_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
hold on;
plot(mean_delay_no_purging,mean_thr_no_purging,'o-c','markerfacecolor','c');
% scatter(mean_delay(1),mean_thr(1),500,'r')
title('Computational load vs Delay - optimum without purging');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');


close(fig6)
fig6 = figure(6); 
scatter(reshape(Dmean_vec(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
hold on;
plot(mean_delay,mean_thr,'o-b','markerfacecolor','b');
% scatter(mean_delay(1),mean_thr(1),500,'r')
title('Computational load vs Delay - optimum with purging');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');

close(fig2)
fig2 = figure(2); 
scatter(reshape(Dmean_vec_gini(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
hold on;
plot(mean_delay_gini,mean_thr_gini,'o-m','markerfacecolor','m');
% scatter(mean_delay(1),mean_thr(1),500,'r')
title('Computational load vs Delay - gini with purging');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');

close(fig1)
fig1 = figure(1); 
scatter(reshape(Dmean_vec_gini_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
hold on;
plot(mean_delay_gini_no_purging,mean_thr_gini_no_purging,'o-g','markerfacecolor','g');
% scatter(mean_delay(1),mean_thr(1),500,'r')
title('Computational load vs Delay - gini without purging');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');

close(fig7)
fig7 = figure(7); 
plot(mean_delay_gini_no_purging,mean_thr_gini_no_purging,'o-g','markerfacecolor','g');
hold on;
plot(mean_delay_gini,mean_thr_gini,'o-m','markerfacecolor','m');
plot(mean_delay_uni_no_purging,mean_thr_uni_no_purging,'o-r','markerfacecolor','r');
plot(mean_delay_uni_with_purging,mean_thr_uni_with_purging,'o-k','markerfacecolor','k');
plot(mean_delay_no_purging,mean_thr_no_purging,'o-c','markerfacecolor','c');
plot(mean_delay,mean_thr,'o-b','markerfacecolor','b');
%scatter(mean_delay(1),mean_thr(1),500,'r')
title('Computational load vs Delay');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');
legend('gini without purging', 'gini with purging','uni without purging', 'uni with purging', 'our without purging', 'our with purging');

close(fig8)
fig8 = figure(8); 
scatter(reshape(Dmean_vec_gini_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
hold on;
scatter(reshape(Dmean_vec_gini(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini(1:niter,:),1,niter*NOmega),10,'m'); %to see the fog
scatter(reshape(Dmean_vec_uni_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_no_purging(1:niter,:),1,niter*NOmega),10,'r'); %to see the fog
scatter(reshape(Dmean_vec_uni_with_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_with_purging(1:niter,:),1,niter*NOmega),10,'k'); %to see the fog
scatter(reshape(Dmean_vec_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_no_purging(1:niter,:),1,niter*NOmega),10,'c'); %to see the fog
scatter(reshape(Dmean_vec(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec(1:niter,:),1,niter*NOmega),10,'b'); %to see the fog
%scatter(mean_delay(1),mean_thr(1),500,'r')
title('Computational load vs Delay');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');
legend('gini without purging', 'gini with purging', 'uni without purging', 'uni with purging', 'our without purging', 'our with purging');
