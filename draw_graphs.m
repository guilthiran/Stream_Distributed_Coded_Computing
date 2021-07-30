function [] = draw_graphs(simu_desc)
%**************************************************************************
%                       UCLOUVAIN/MIT CONFIDENTIAL
%                           ALL RIGHTS RESERVED
%**************************************************************************

%**************************************************************************
% Generate plot for the decentralized computing scheme
%
%
% Author:           Guillaume Thiran    (guillaume.thiran@uclouvain.be)
%                   2021/01/08
% Last modified :   Guillaume Thiran    (guillaume.thiran@uclouvain.be)
%                   2020/01/31
%**************************************************************************
    
%1) graph of the total time
figure(1); 
plot(simu_desc.analysis.total_time_job); 
xlabel('Job index')
ylabel('(Service + Waiting) time [slots]');
title('Sojourn time of the jobs');

%2) graph of estimation
figure(2);
%estimation
subplot(121)
tmax = max(simu_desc.results.index_service);
for p=1:simu_desc.parameters.P
    plot(linspace(1,tmax,simu_desc.results.index_service(p)),simu_desc.results.m1_est(p,1:simu_desc.results.index_service(p))','b');
    hold on;
end
%true value
plot(repmat([1;tmax],1,simu_desc.parameters.P),[simu_desc.results.m1_est(:,1)';simu_desc.results.m1_est(:,1)'],'k')

subplot(122)
for p=1:simu_desc.parameters.P
    plot(linspace(1,tmax,simu_desc.results.index_service(p)),simu_desc.results.m2_est(p,1:simu_desc.results.index_service(p))','b');
    hold on;
end
%true value
plot(repmat([1;tmax],1,simu_desc.parameters.P),[simu_desc.results.m2_est(:,1)';simu_desc.results.m2_est(:,1)'],'k')

%3) graph of the worker's use
% n_bins =10; 
% bins = floor(linspace(1, simu_desc.results.index_task-1,n_bins+1));
% w_vec = zeros(simu_desc.parameters.P,n_bins+1); 
% 
% for b=1:n_bins
%     for m=bins(b):bins(b+1)
%         p=simu_desc.results.task_desc(m).worker;
%         w_vec(p,b)=w_vec(p,b)+1;
%     end
% end
% [phi,~] = optimum_split(simu_desc.parameters.r_comp_scaled,simu_desc.parameters.r_comm_scaled,simu_desc.parameters.phi_l);
% w_vec(:,end)=phi';
% w_vec=w_vec./repmat(sum(w_vec,1),simu_desc.parameters.P,1);
% figure(3);
% bar([1:n_bins n_bins+2],w_vec','stacked')
% xlabel('Time');
% ylabel('Worker use');
% title('Worker use as a function of time(right bar is the optimum one)');
% 
% %graph of master queue length
% figure(4)
% plot(simu_desc.results.length_queue(1:simu_desc.results.t),'b','linewidth',1.5);
% xlabel('time');
% title('Master queue length');
% 
%graph of workers' queue length
figure(5)
subplot(131)
plot(simu_desc.results.Queue_length_in(:,1:simu_desc.results.t)','linewidth',1.5);
xlabel('time');
title('Comm in queue length');
subplot(132)
plot(simu_desc.results.Queue_length_comp(:,1:simu_desc.results.t)','linewidth',1.5);
xlabel('time');
title('Comp queue length');
subplot(133)
plot(simu_desc.results.Queue_length_out(:,1:simu_desc.results.t)','linewidth',1.5);
xlabel('time');
title('Comm out queue length');

end