clearvars;
close all;

%load('Distributed_computing_trdoff_simu22021_3_1_5_48');
load('simu_results/Distributed_computing_data_non_iterative_paper');

close(fig3)
fig3 = figure(3); 

for o = 1:NOmega
    idx_neg = Dmean_vec_uni_no_purging(:,o)<=350;
    mean_replace = mean(Dmean_vec_uni_no_purging(~idx_neg,o));
    Dmean_vec_uni_no_purging(idx_neg,o)=mean_replace;
    
    idx_leq1 = glo_comp_load_ratio_vec_uni_no_purging(:,o)<1;
    repl = mean(glo_comp_load_ratio_vec_uni_no_purging(~idx_leq1,o));
    glo_comp_load_ratio_vec_uni_no_purging(idx_leq1,o)=repl;
end

mean_delay_uni_no_purging = mean(Dmean_vec_uni_no_purging(1:niter,:),1);
mean_thr_uni_no_purging = mean(glo_comp_load_ratio_vec_uni_no_purging(1:niter,:),1);
scatter(reshape(Dmean_vec_uni_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
hold on;
plot(mean_delay_uni_no_purging,mean_thr_uni_no_purging,'o-r','markerfacecolor','r');
%scatter(mean_delay_uni_no_purging(1),mean_thr_uni_no_purging(1),500,'r')
title('Computational load vs Delay - uniform without purging');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');


for o = 1:NOmega
    idx_neg = Dmean_vec_uni_with_purging(:,o)<=350;
    mean_replace = mean(Dmean_vec_uni_with_purging(~idx_neg,o));
    Dmean_vec_uni_with_purging(idx_neg,o)=mean_replace;
    
    idx_leq1 = glo_comp_load_ratio_vec_uni_with_purging(:,o)<1;
    repl = mean(glo_comp_load_ratio_vec_uni_with_purging(~idx_leq1,o));
    glo_comp_load_ratio_vec_uni_with_purging(idx_leq1,o)=repl;
end

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


for o = 1:NOmega
    idx_neg = Dmean_vec_no_purging(:,o)<=350;
    mean_replace = mean(Dmean_vec_no_purging(~idx_neg,o));
    Dmean_vec_no_purging(idx_neg,o)=mean_replace;
    
    idx_leq1 = glo_comp_load_ratio_vec_no_purging(:,o)<1;
    repl = mean(glo_comp_load_ratio_vec_no_purging(~idx_leq1,o));
    glo_comp_load_ratio_vec_no_purging(idx_leq1,o)=repl;
end


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

for o = 1:NOmega
    idx_neg = Dmean_vec(:,o)<=350;
    mean_replace = mean(Dmean_vec(~idx_neg,o));
    Dmean_vec(idx_neg,o)=mean_replace;
    
    idx_leq1 = glo_comp_load_ratio_vec(:,o)<1;
    repl = mean(glo_comp_load_ratio_vec(~idx_leq1,o));
    glo_comp_load_ratio_vec(idx_leq1,o)=repl;
end

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

for o = 1:NOmega
    idx_neg = Dmean_vec_gini(:,o)<=350;
    mean_replace = mean(Dmean_vec_gini(~idx_neg,o));
    Dmean_vec_gini(idx_neg,o)=mean_replace;
    
    idx_leq1 = glo_comp_load_ratio_vec_gini(:,o)<1;
    repl = mean(glo_comp_load_ratio_vec_gini(~idx_leq1,o));
    glo_comp_load_ratio_vec_gini(idx_leq1,o)=repl;
end

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

for o = 1:NOmega
    idx_neg = Dmean_vec_gini_no_purging(:,o)<=350;
    mean_replace = mean(Dmean_vec_gini_no_purging(~idx_neg,o));
    Dmean_vec_gini_no_purging(idx_neg,o)=mean_replace;
    
    idx_leq1 = glo_comp_load_ratio_vec_gini_no_purging(:,o)<1;
    repl = mean(glo_comp_load_ratio_vec_gini_no_purging(~idx_leq1,o));
    glo_comp_load_ratio_vec_gini_no_purging(idx_leq1,o)=repl;
end

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

fig9=figure(9);
%semilogy(Omega_vec,mean_delay_gini_no_purging,'o-g','markerfacecolor','g');
semilogy(Omega_vec,mean_delay_gini,'o-m','markerfacecolor','m');
hold on;
%semilogy(Omega_vec,mean_delay_uni_no_purging,'o-r','markerfacecolor','r');
semilogy(Omega_vec,mean_delay_uni_with_purging,'o-k','markerfacecolor','k');
%semilogy(Omega_vec,mean_delay_no_purging,'o-c','markerfacecolor','c');
semilogy(Omega_vec,mean_delay,'o-b','markerfacecolor','b');
%scatter(mean_delay(1),mean_thr(1),500,'r')
title('Delay vs Omega');
ylabel('Mean delay');
xlabel('Omega');
legend('gini with purging', 'uni with purging', 'our with purging');
%legend('gini without purging', 'gini with purging','uni without purging', 'uni with purging', 'our without purging', 'our with purging');


close(fig8)
fig8 = figure(8); 
%scatter(reshape(Dmean_vec_gini_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
%scatter(reshape(Dmean_vec_gini(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini(1:niter,:),1,niter*NOmega),10,'m'); %to see the fog
scatter(reshape(Dmean_vec_uni_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_no_purging(1:niter,:),1,niter*NOmega),10,'r'); %to see the fog
hold on;
scatter(reshape(Dmean_vec_uni_with_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_uni_with_purging(1:niter,:),1,niter*NOmega),10,'k'); %to see the fog
scatter(reshape(Dmean_vec_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_no_purging(1:niter,:),1,niter*NOmega),10,'c'); %to see the fog
scatter(reshape(Dmean_vec(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec(1:niter,:),1,niter*NOmega),10,'b'); %to see the fog
%scatter(mean_delay(1),mean_thr(1),500,'r')
title('Computational load vs Delay');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');
legend('uni without purging', 'uni with purging', 'our without purging', 'our with purging');

close(fig8)
fig8 = figure(8); 
nn=100;
bound = 2000;
D1 = reshape(Dmean_vec_uni_no_purging(1:nn,:),1,nn*NOmega);
G1 = reshape(glo_comp_load_ratio_vec_uni_no_purging(1:nn,:),1,nn*NOmega);
idxok = D1<bound;
scatter(D1(idxok),G1(idxok),10,'r'); %to see the fog

%scatter(reshape(Dmean_vec_gini_no_purging(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini_no_purging(1:niter,:),1,niter*NOmega),10,'g'); %to see the fog
%scatter(reshape(Dmean_vec_gini(1:niter,:),1,niter*NOmega),reshape(glo_comp_load_ratio_vec_gini(1:niter,:),1,niter*NOmega),10,'m'); %to see the fog
hold on;
D1 = reshape(Dmean_vec_uni_with_purging(1:nn,:),1,nn*NOmega);
G1 = reshape(glo_comp_load_ratio_vec_uni_with_purging(1:nn,:),1,nn*NOmega);
idxok = D1<bound;
scatter(D1(idxok),G1(idxok),10,'k'); %to see the fog
D1 = reshape(Dmean_vec_no_purging(1:nn,:),1,nn*NOmega);
G1 = reshape(glo_comp_load_ratio_vec_no_purging(1:nn,:),1,nn*NOmega);
idxok = D1<bound;
scatter(D1(idxok),G1(idxok),10,'c'); %to see the fog
D1 = reshape(Dmean_vec(1:nn,:),1,nn*NOmega);
G1 = reshape(glo_comp_load_ratio_vec(1:nn,:),1,nn*NOmega);
idxok = D1<bound;
scatter(D1(idxok),G1(idxok),10,'b'); %to see the fog
%scatter(reshape((1:niter,:),1,niter*NOmega),reshape((1:niter,:),1,niter*NOmega),10,'b'); %to see the fog
%scatter(mean_delay(1),mean_thr(1),500,'r')
title('Computational load vs Delay');
xlabel('Mean delay');
ylabel('Computational load (number of tasks/optimum number)');
legend('uni without purging', 'uni with purging', 'our without purging', 'our with purging');


drawnow
