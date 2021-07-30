function simu_desc_upd = analysis_results(simu_desc, verbose)
%**************************************************************************
%                       UCLOUVAIN/MIT CONFIDENTIAL
%                           ALL RIGHTS RESERVED
%**************************************************************************

%**************************************************************************
% Analysis of the results if the decentralized streaming computation scheme
% This function prints information if verbose=1 and updates the simu_desc
% structure
%
%
% Author:           Guillaume Thiran    (guillaume.thiran@uclouvain.be)
%                   2021/01/08
% Last modified :   Guillaume Thiran    (guillaume.thiran@uclouvain.be)
%                   2020/02/01
%**************************************************************************

%1) in order time spent to solve jobs
total_time_job = zeros(1,simu_desc.parameters.J); %starting at the moment the job arrives
last_job = 0; 
for j=1:simu_desc.parameters.J
    end_timej = simu_desc.results.job_desc(j).end_time;
    if end_timej<=last_job
        end_timej=last_job;
    else
        last_job=end_timej;
    end
    total_time_job(j) = end_timej-simu_desc.results.job_desc(j).arrival_time;
end
mean_total_time = mean(total_time_job+simu_desc.parameters.Denc+simu_desc.parameters.Ddec);
max_total_time = max(total_time_job+simu_desc.parameters.Denc+simu_desc.parameters.Ddec);
%optimum time if queues at the workers
%[phi,overload] = optimum_split(simu_desc.parameters.r_comp_scaled,simu_desc.parameters.r_comm_scaled,simu_desc.parameters.phi_l);
% if overload %if even with queues the arrival rate is too high
%     warning('Arrival rate of jobs exceed optimum global service rate');
% end
%the optimum total time is computed as the mean total time per task, 
%times the total number of tasks per jobs.   
% opt_total_time = simu_desc.parameters.Denc+simu_desc.parameters.Ddec+1/simu_desc.parameters.P*simu_desc.parameters.K*simu_desc.parameters.Omega*sum(phi./(simu_desc.parameters.mu_rate_vec_task'-phi*simu_desc.parameters.lambda_J_task))+1/simu_desc.parameters.P*sum(phi*simu_desc.parameters.I./simu_desc.parameters.c_rate_vec');
% if verbose
%     fprintf('waiting + service time. mean: %d; max: %3.1d [slots]\n',round(mean_total_time) ,max_total_time);
%     fprintf('Optimum mean total time %d [slots] \n',round(opt_total_time));
% end

%2) global service rate
% mu_tot= 1/mean(active_time_job); %tasks divided by time
% if verbose
%     fprintf('Global service rate %1.5f (arrival rate %1.5f) [job/slots] \n',mu_tot ,simu_desc.parameters.lambda_J);
% end
% if mu_tot < simu_desc.parameters.lambda_J %if the global service rate is too low, the results will be useless.
%                      %As it will be impossible to define a "service time"
%                      %and as when more jobs will be added, the total time
%                      %will just go to infinity, and will not converge to a
%                      %mean value. 
%     warning('Global service rate lower than jobs arrival rate. Results may be unreliable.');
% end

%3) Computationnal load. Due to the feedback delay, more tasks will be
%computed than what we need. Hence, it is important to monitor the
%computationnal load of the system. For this, we need to consider tasks
%that either succeed, or fail due to time-out.
comp_load_vec=zeros(1,simu_desc.parameters.J);
for j=1:simu_desc.parameters.J %for all tasks that have effectively been assigned
    comp_load_vec(j) = simu_desc.results.job_desc(j).n_OK_tasks;
end
%comp_load_vec_ratio = comp_load_vec./simu_desc.parameters.K; %computational load with regards to the optimal one
glo_comp_load = sum(comp_load_vec); %global computational load
glo_comp_load_ratio=glo_comp_load./(simu_desc.parameters.J*simu_desc.parameters.K); %global ratio
if verbose
    fprintf('Computational load ratio %2.3f \n',glo_comp_load_ratio);
end

%4) Communication load. We take into account communications between master
%and workers, and between workers and fusion node.
% comm_load_vec=zeros(1,simu_desc.parameters.J);
% for m=1:simu_desc.results.index_task_worker %for all tasks that have effectively been assigned
%     j=simu_desc.results.task_desc(m).job;
%     %1 comm use if fail before or fail after
%     %2 if successfull task
%     comm_load_vec(j) = comm_load_vec(j) +simu_desc.results.task_desc(m).fail_before+simu_desc.results.task_desc(m).fail_after+2*simu_desc.results.task_desc(m).sucess;
% end
% comm_load_vec_ratio = comm_load_vec./(2*simu_desc.parameters.K); %communication load with regards to the optimal one
%The optimal one is given by two times the recovery threshold as if there
%are no failues, this means we need K communications between master and
%workers and between workers and fusion.
% glo_comm_load = sum(comm_load_vec); %global communication load
% glo_comm_load_ratio=glo_comm_load./(simu_desc.parameters.J*2*simu_desc.parameters.K); %global ratio
% if verbose
%     fprintf('Communication load ratio %2.3f \n',glo_comm_load_ratio);
% end

%5) worker use. As we will pay for the use of the workers, we do not want
%them to be idle. Hence, we monitor the percentage of time during which
%they are idle.
% idle_per_worker = sum(simu_desc.results.current_task(:,1:(simu_desc.results.t-1))==0,2);
% percentage_per_worker = idle_per_worker/simu_desc.results.t*100;
% percentage_idle_tot = sum(idle_per_worker)/(simu_desc.results.t*simu_desc.parameters.P)*100;
% if verbose
%     fprintf('Idle percentage %2.3f %% \n',percentage_idle_tot);
% end

%6) For simulation purposes, analysis of the choice of max task/max time
if verbose
    fprintf('Task index %d / %d \n', simu_desc.results.index_task, simu_desc.parameters.M_init);
    fprintf('Time index %d / %d \n', simu_desc.results.t, simu_desc.parameters.T);
end

% %7) Throughput: i.e. number of tasks divided by optimal one. 
% opt_n_tasks = simu_desc.parameters.J*simu_desc.parameters.K;
% n_tasks = simu_desc.results.index_task;
% ratio_throughput = n_tasks/opt_n_tasks;
% if verbose
%     fprintf('Throughput ratio %2.3f \n', ratio_throughput);
% end
    
%------------------------results' analysis------------------------%
%simu_desc.analysis.active_time_job=active_time_job; %time between beginning of processing and end
simu_desc.analysis.total_time_job=total_time_job; %total time between arrival and end
simu_desc.analysis.mean_total_time = mean_total_time; %mean total time
simu_desc.analysis.max_total_time = max_total_time; %max total time
simu_desc.analysis.glo_comp_load_ratio=glo_comp_load_ratio;
%simu_desc.analysis.mu_tot=mu_tot; %whole system's service rate
%simu_desc.analysis.phi_opt=phi; %optimum split if workers have queues
%simu_desc.analysis.opt_total_time=opt_total_time; %optimum mean response time if workers have queues

% simu_desc.analysis.comp_load_vec=comp_load_vec; %computational load
% simu_desc.analysis.comp_load_vec_ratio=comp_load_vec_ratio; %ratio of computational load
% simu_desc.analysis.glo_comp_load = glo_comp_load; %global computational load
% simu_desc.analysis.glo_comp_load_ratio=glo_comp_load_ratio; %global computational load ratio
% 
% simu_desc.analysis.comm_load_vec=comm_load_vec; %communication load
% simu_desc.analysis.comm_load_vec_ratio=comm_load_vec_ratio; %ratio of communication load
% simu_desc.analysis.glo_comm_load = glo_comm_load; %global computational load
% simu_desc.analysis.glo_comm_load_ratio=glo_comm_load_ratio; %global computational load ratio
% 
% simu_desc.analysis.percentage_per_worker=percentage_per_worker; %idle percentage per worker
% simu_desc.analysis.percentage_idle_tot=percentage_idle_tot; %idle percentage
% 
% simu_desc.analysis.ratio_throughput=ratio_throughput; %throughput ratio

simu_desc_upd=simu_desc;
end

