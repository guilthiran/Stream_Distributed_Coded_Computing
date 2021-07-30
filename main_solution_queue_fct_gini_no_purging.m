function simu_desc = main_solution_queue_fct_gini_no_purging(parameters,saved,draw,verbose)

%**************************************************************************
%                       UCLOUVAIN/MIT CONFIDENTIAL
%                           ALL RIGHTS RESERVED
%**************************************************************************

%**************************************************************************
% Streaming distributed computing
%
% Simulation of a decentralized computing scheme 
% Workers are equipped with queues
%
% Author:           Guillaume Thiran    (guillaume.thiran@uclouvain.be)
%                   2021/12/15
% Last modified :   Guillaume Thiran    (guillaume.thiran@uclouvain.be)
%                   2021/02/01
%**************************************************************************


%% Parameters
simu_desc=struct; 
parameters.verbose=verbose; %1 if information is printed
parameters.save=saved; %1 if information is saved
parameters.draw=draw; %1 if plots are generated

%------------------------Initial choice of workers-------------------%
parameters = no_choose_workers(parameters);

%------------------------System model (based on input parameters)------------------------%
parameters.interarrival_time=exprnd(1/parameters.lambda_J,1,parameters.J-1); %interrival time of the jobs, exponentially distributed
parameters.arrival_time_vec = [1, 1+cumsum(parameters.interarrival_time)]; %arrival time in slots (with first job arriving at t=1)

parameters.fact_task=3; %security factor for the tasks (in terms of number of task we expect)
parameters.M=ceil(parameters.fact_task * parameters.K*parameters.J); %total number of tasks
M_init=parameters.M;
parameters.M_init = M_init; %initial estimation of the number of tasks
parameters.fact_time=2; %security factor for time
parameters.T=ceil(parameters.fact_time*parameters.M/sum(parameters.mu_rate_vec_task)); %maximum time
parameters.tt=1:parameters.T;

% parameters.t_feedback=[0 parameters.tt]; %feedback time slots, denoting the time slots where a feedback is sent 
% parameters.t_feedback=[parameters.t_feedback parameters.t_feedback(end)+parameters.feedback_interval]; %for simulation, we need this additionnal feedback time (which will not be useful as we will reach the max time before)
% parameters.feedback_delay=1; %feedback delay (from fusion to master), in [slots] 


%shifted exponential
parameters.service_time=zeros(parameters.P,parameters.T); %service time of worker p is task assigned at time slot t
for p=1:parameters.P
    parameters.service_time(p,:)=exprnd(1/parameters.mu_rate_vec_task(p),1,parameters.T);
end
% parameters.mu_est_init=parameters.mu_rate_vec_task;%initial estimation
% lambda_J_task = parameters.lambda_J_task; %arrival rate in task/time slot

% parameters.t_out_vec=(1+parameters.fact_time_out)./parameters.mu_rate_vec_task; 

simu_desc.parameters=parameters;

%% Initialization (Normally, no parameter is present from here to the bottom of the file)

%------------------------Jobs------------------------%
J = parameters.J; %Number of jobs
arrival_time_vec = parameters.arrival_time_vec;%arrival time in slots (with first job arriving at t=1)
%create a structure containing a description of the jobs
job_desc = struct;
for j=1:J
    job_desc(j).rec_thres = parameters.K; 
    job_desc(j).arrival_time = ceil(arrival_time_vec(j)); %arrival time ceiled to the following time slot 
%     job_desc(j).n_sent_tasks = 0; %number of tasks that have been sent to a worker for the job
%     job_desc(j).n_sent_tasks_prev=0;
%     job_desc(j).fail=0;
    job_desc(j).n_OK_tasks = 0; %number of tasks that have suceeded (as observed at the fusion node)
%     job_desc(j).n_OK_tasks_master = 0;
    job_desc(j).state = 0; %0 when not available, 
                           %1 when available but not started, 
                           %2 when started but not finished, 
                           %3 when finished at the fusion node 
                           %4 when finished and known at the master node  
    job_desc(j).end_time = 0; %time at which a job is finished at the fusion node 
%     job_desc(j).n_OK_fail=0;
end
jobs_all_received=0; %1 when all jobs are received by the master, waiting to be processed

%------------------------Tasks------------------------%
M = parameters.M; %maximum number of tasks, obtained as the sum of the recovery thres times a factor
%create a structure containing a description of the tasks
task_desc = struct;
%t_out_vec=parameters.t_out_vec;
for m=1:M
    task_desc(m).job = 0; %job which corresponds to the task
%     task_desc(m).worker = 0; %worker processing the task
%     task_desc(m).sucess = 0; %1 if the task has succeeded 
    task_desc(m).t_success = 0; %if sucess, time of sucess, in slot
%     task_desc(m).state_master=0;
    task_desc(m).t_tasks=0;
end

%------------------------Workers------------------------%
P=parameters.P; %number of workers
D_comm_in_vec_task = parameters.D_comm_in_vec_task; %input communication delay
D_comm_out_vec_task = parameters.D_comm_out_vec_task; %output communication delay
 %input communication delay
% t_out_vec=parameters.t_out_vec;

%------------------------Time------------------------%
T = parameters.T; %the maximum number of time slots is computed as the maximum number of tasks divided by the total service rate, with a security factor
%if the arrival rate of the job is very low, it may happen the above
%maximum time is not enough. 
% t_feedback =  parameters.t_feedback;%time slots for which a feedback is sent (not including the feedback delay)
%                                                           %When a feedback is sent at time t, we assume we know everything up to time t, included 
% feedback_delay=parameters.feedback_delay; %delay between a feedback is sent and the feedback is received
                                                        
%------------------------Workers (bis)------------------------%
service_time = parameters.service_time; 

Queue_main = [];
Queue_main_length=0;

Queue_length_in = zeros(P,T);
Queue_comm_in =  cell(P,1);
comm_in = -1*ones(P,2);
Queue_length_comp = zeros(P,T);
Queue_comp = cell(P,1);
comp = -1*ones(P,3);
Queue_length_out = zeros(P,T);
Queue_comm_out =  cell(P,1);
comm_out = -1*ones(P,2);

%------------------------Master------------------------%
job_finished = zeros(1,J); %1 if job j is known to be finished
max_time =0; %index of the time at which we known the last job is finished
% r_comm_scaled = parameters.r_comm_scaled;
%r_comp_scaled = parameters.r_comp_scaled;
% phi_l=parameters.phi_l;
Omega_forward=parameters.Omega_forward;
% index_service = zeros(P,1); %up to which index service time have been collected
% service_time_collect=zeros(P,T); %collection of the service time
% mu_est=zeros(P,length(t_feedback));%estimation of the service rate of the workers at each time step
% mu_est(:,1)=parameters.mu_est_init; %initialization
% beta=parameters.beta; %memory factor for update of mu
% job_finished_master=-1*ones(1,J);
% tprime = zeros(J,1); %sent time of the last task we have received feedback for
% 
% jminj = 1; %index of the first (oldest) job still being processed
% jmaxj = 0; %index of the last (newest) job being processed.
% P_failure = mean(exp(- parameters.mu_rate_vec_task.*t_out_vec));
% m_min = 0; %oldest task for which we know for sure either it
%                                     %has succeded or not

%-------------------Fusion---------------------%
task_received = zeros(P,T); %The element in (p,t) gives the ID of the task which is finished at worker p for time t, if any


%% distributed computation simulation
%to know at which index looking in the different vectors
index_t_job_arrival = 1;
index_task=1; %index of the next task to be scheduled
t=1; %time step
end_simu=0; %1 if all jobs are finished

last_feedback=0;

index_t_feedback = 2; %the vector starts at 0 to simplify the algorithmic part (we are looking for the previous feedback, also avoid index problems at the end)
while ~end_simu && t<=T %at each time slot, if we do not exceed the maximum number of task
%      feedback_arrived=0;
%      %feedback
%      %------------------------Feedback------------------------%
%     tf_sent = t_feedback(index_t_feedback); %time slot of the next feedback (when sent)
%     tf_rec = tf_sent+feedback_delay; %time slot of the next feedback (when received)
%     if t==tf_rec %if feedback time 
%         feedback_arrived=1; 
%         
%         %We need to look to tasks received between t_feedback(index_t_feedback-1)+1 and tf_sent
%         t_begin = t_feedback(index_t_feedback-1)+1;
%         new_feedbacks = task_received(:,t_begin:tf_sent); %tasks which are concerned by this feedback
%         
%         %Use task information to update job and parameters estimation
%         %By default, same number of received task
%         completed_tasks = sort(nonzeros(new_feedbacks)); 
%         worker_last=zeros(1,P); %index of the last task succeeded by the worker
%         if ~isempty(completed_tasks)
%             for m=completed_tasks' %for each completed tasks
%                 task_desc(m).state_master=1;
%                 %Update the jobs' informations (i.e. the number of tasks the master know are received,
%                 %and the number of missing degree of freedom)
%                 j = task_desc(m).job;
%                 job_desc(j).n_OK_tasks_master=job_desc(j).n_OK_tasks_master+1;
%                 tprime(j)=max(tprime(j),task_desc(m).t_assignment);
%                 
%                 %uses the time of arrival to update the estimated rate of the workers
%                 p=task_desc(m).worker;
%                 index_service(p)=index_service(p)+1;
%                 t_serv = task_desc(m).t_tasks;%task_desc(m).t_success-task_desc(m).t_assignment;
%                 service_time_collect(p,index_service(p))=t_serv; %collection of the service time
%                 
%                 %update the lastt ask the worker has processed
%                 worker_last(p)=m;
%             end
%             
%             %check job
%             jminnew=jminj;
%             for j=jminj:jmaxj
%                 if job_desc(j).n_OK_tasks_master>=job_desc(j).rec_thres
%                     job_finished_master(j)=1;
%                     if j==jminnew
%                         jminnew=j+1;
%                     end
%                 end
%             end
%             jminj=jminnew;
% 
%             %At this point, we have not received the feedback for some of
%             %the tasks. We need to determine which tasks are finished,
%             %based on the fact that we have
%             %received a feedback for a task the was assigned to same worker
%             %after
%             m = m_min+1; %m_min+1 denotes the oldest take for which we are not sure yet
%             m_min_set=0;
%             while m<=index_task-1 && task_desc(m).t_assignment <tf_sent
%                 p=task_desc(m).worker; %worker of the task
%                 if task_desc(m).state_master==0 %if unknown state
%                     %if a more recent task is finished for this worker or
%                     %if the interval of time between the assignment of the
%                     %task and the time of feedback transmission is greater
%                     %than the time-out ==> task has failed
%                     if m < worker_last(p)
%                         %task has failed
%                         %task_desc(m).state_master=2;
%                         %update info
%                         j=task_desc(m).job;
%                         job_desc(j).fail = job_desc(j).fail+1;
%                         task_desc(m).state_master=1;
%                     elseif ~m_min_set %if we do not know the state and m_min has not been updated yet
%                         m_min=m-1;
%                         m_min_set=1;
%                     end
%                 end
%                 m=m+1;
%             end
%             if ~m_min_set %if m_min has not been updated, this means we now the state of all tasks
%                 m_min=m-1; %m_min equal the index of the last task for which we have received feedback
%             end
% 
%             %change mu_est: for the moment, censored estimation
%             for p=1:P
%                 r=index_service(p);
%                 x = sum(service_time_collect(p,1:index_service(p)));
%                 if r>0
%                     x0=r/x; %brute ML estimator
%                     mu_est(p,index_t_feedback)=fzero(@(la) trunc_esti_exp(la,x,r,t_out_vec(p)),x0);
%                     if ~(mu_est(p,index_t_feedback)>0) %at the beginning, this may fail.
%                         mu_est(p,index_t_feedback)=r/x; %in this case, brute ML
%                     end
%                 else
%                     mu_est(p,index_t_feedback)=0; %if only time-out  
%                 end
%                 %use memory in order to avoid bad estimation 
%                  mu_est(:,index_t_feedback) = mu_est(:,index_t_feedback-1)*(1-beta)+beta*mu_est(:,index_t_feedback);
% 
%             end
%         else
%             %use memory in order to avoid bad estimation 
%             mu_est(:,index_t_feedback) = mu_est(:,index_t_feedback-1);
% 
%         end
%         
%         index_t_feedback=index_t_feedback+1; %index for the next feedback time
%     end
    
    
    
    %Reverse direction. Start with Comm out
    %------------------------Comm_out------------------------%
    jmin = J; %to remember which job should be checked
    jmax = 0;
    for p=1:P
        if comm_out(p,1)>=0 && comm_out(p,2)<=t %if task in comm, and comm is finished
            m=comm_out(p,1); %task index
            j=task_desc(m).job;
            job_desc(j).n_OK_tasks = job_desc(j).n_OK_tasks+1;
            
            task_received(p,t)=m; 
            
            jmin = min(jmin,j);
            jmax = max(jmax,j);
            
            comm_out(p,:)=-1;
        end
    end
    %for all jobs, check whether they are finished
    for j = jmin:jmax
        if job_desc(j).state<4 && job_desc(j).n_OK_tasks>=job_desc(j).rec_thres %if job is not already finished and if it is now finished
            job_desc(j).state=4; %job is finished 
            job_desc(j).end_time = t;
            job_finished(j)=1;
            max_time=t;
            
          %   [comm_in,Queue_comm_in,comp,Queue_comp,comm_out,Queue_comm_out]...
          %=purge_queues(comm_in,Queue_comm_in,comp,Queue_comp,comm_out,Queue_comm_out,j,P,task_desc);
            
        end
    end
     %     
    %------------------------Queue comm out------------------------%
    for p=1:P
        if comm_out(p,1)<0 %if no task currently transmitted
            if ~isempty(Queue_comm_out{p})
                comm_out(p,1)= Queue_comm_out{p}(1);
                comm_out(p,2) = t+D_comm_out_vec_task(p);
                
                Queue_comm_out{p}(1)=[];
            end
        end
        Queue_length_out(p,t)=length(Queue_comm_out{p});
    end
    
    %------------------------Computation--------------------------%
    for p=1:P
        if comp(p,1)>=0 && comp(p,2)<=t %if task in comp, and comp is finished
            if comp(p,3) == -1 %normal task
                m=comp(p,1); %task index
                Queue_comm_out{p}=[Queue_comm_out{p} m];
            end 
            comp(p,:)=-1;
        end
    end
    
    %------------------------Queue comp------------------------%
    for p=1:P
        if comp(p,1)<0 %if no task currently computed
            if ~isempty(Queue_comp{p})
                m = Queue_comp{p}(1);
                comp(p,1) =m;
                t_tasks = service_time(p,t);
                task_desc(m).t_tasks=t_tasks;
                comp(p,2) = t+t_tasks;
                
                Queue_comp{p}(1)=[];
            end
        end
        Queue_length_comp(p,t)=length(Queue_comp{p});
    end
    
    %------------------------Communication in--------------------------%
    for p=1:P
        if comm_in(p,1)>=0 && comm_in(p,2)<=t %if task in comm in, and comm in is finished
            m=comm_in(p,1); %task index
            Queue_comp{p}=[Queue_comp{p} m];
            
            comm_in(p,:)=-1;
        end
    end
    
    %------------------------Queue comm in------------------------%
    for p=1:P
        if comm_in(p,1)<0 %if no task currently transmitted
            if ~isempty(Queue_comm_in{p})
                comm_in(p,1)= Queue_comm_in{p}(1);
                comm_in(p,2) = t+D_comm_in_vec_task(p);
                
                Queue_comm_in{p}(1)=[];
            end
        end
        Queue_length_in(p,t)=length(Queue_comm_in{p});
    end   
    
    %---------------non-causal assignment---------------------------%
    T_ass = floor(comp(:,2)-D_comm_in_vec_task');
    
    for p=1:P   
        if (isempty(Queue_comm_in{p}) && comm_in(p,2)<0 && isempty(Queue_comp{p}) && comp(p,2)<0)...
                ||t==T_ass(p) %if worker will be idle after comm 
            if Queue_main_length>0
                m = Queue_main(1);
                Queue_comm_in{p} = [Queue_comm_in{p} m];
                Queue_main_length=Queue_main_length-1;
                Queue_main(1)=[];
            end
        end
    end
    
    %-----------------------Supportive_tasks---------------------------%
%     nJobs = jmaxj-jminj+1;
%     Delta_mat = zeros(1,nJobs);
%     if feedback_arrived
%         %compute the second Delta for each job
%         idxj=1;
%         for j = jminj:jmaxj
% 
%             interval=t-last_feedback;
%             n_OK_fail = job_desc(j).n_OK_fail;
%             n_fails = job_desc(j).fail;
%             Delta_mat(1,idxj) = max(n_fails-n_OK_fail,0)/interval;%max(( min(num*(1-P_failure),K)-Bt)/interval,0);
% 
%             job_desc(j).n_sent_tasks_prev=job_desc(j).n_sent_tasks;
% 
%             
%             idxj=idxj+1;
%         end
%         last_feedback=t;
%             
%     end
    
    %------------------------Job arrival & forwarding------------------------%
    %while loop as it is possible two jobs arrive at the same time
    %optimum_split
%     r_comp_scaled_est = mu_est(:,index_t_feedback-1)'/lambda_J_task;
%     [phi,~] = optimum_split(r_comp_scaled_est,r_comm_scaled,phi_l);
%         
    while ~jobs_all_received && t==job_desc(index_t_job_arrival).arrival_time
        %forward $phi_pK\OmegaTilde tasks to each worker$
        j=index_t_job_arrival;
%         job_finished_master(j)=0;
%         jmaxj=max(jmaxj,j);
        n_tasktot = ceil(job_desc(j).rec_thres*Omega_forward);
        list_ID_tasks = index_task:index_task+n_tasktot-1; %ID of the tasks
        index_task=index_task+n_tasktot;
        Queue_main = [Queue_main list_ID_tasks]; %#ok<AGROW>
        Queue_main_length=Queue_main_length+n_tasktot;
        
        if index_task>M
            for m=(M+1):(2*M)
                task_desc(m).job = 0; %job which corresponds to the task
                %     task_desc(m).worker = 0; %worker processing the task
                %     task_desc(m).sucess = 0; %1 if the task has succeeded 
                    task_desc(m).t_success = 0; %if sucess, time of sucess, in slot
                %     task_desc(m).state_master=0;
                    task_desc(m).t_tasks=0;  
            end
            M=2*M; %update the maximum number of tasks
            simu_desc.parameters.M=M; %update in the simu desc
        end

        for m=list_ID_tasks %task index
            %update task information
            task_desc(m).job = j; %job which corresponds to the task
            %task_desc(m).t_assignment = t; %assignment time, in slot
            %task_desc(m).worker = p; %worker processing the task
        end
        
        %if all jobs are done, we do not need to check the arrival
        if index_t_job_arrival==J
            jobs_all_received=1;
        end
        index_t_job_arrival=index_t_job_arrival+1;
    end
        
%         n_tasks = ceil(phi*job_desc(j).rec_thres*Omega_forward); %number of tasks per worker
        
%         job_desc(j).n_sent_tasks=job_desc(j).n_sent_tasks+sum(n_tasks);
%         job_desc(j).n_OK_fail  = sum(n_tasks)-job_desc(j).rec_thres;
%             
    
  
    
    
   
    
    %------------------------Stopping criterion------------------------%
    %If all jobs are known to be finished, stop
    if sum(job_finished)==J
        end_simu=1; 
    end
    t=t+1;
end

%------------------------Simulation results------------------------%
results.job_desc=job_desc; %description of the jobs
results.job_finished=job_finished; %1 if job j is finished
results.task_desc=task_desc; %description of the tasks
results.t=t; %time index
results.max_time = max_time; %time at which the last job has been received by the fusion node
results.index_task=index_task; %task index
results.index_t_job_arrival=index_t_job_arrival; %job arrival index
results.end_simu=end_simu; %1 if the simu has correctly finished, 0 if maximum time is reached
results.Queue_length_in=Queue_length_in;
results.Queue_length_comp=Queue_length_comp;
results.Queue_length_out=Queue_length_out;
simu_desc.results=results;

%% Sanity check
%sanity_check(simu_desc,parameters.verbose);

%% Analysis of the results
simu_desc = analysis_results(simu_desc, parameters.verbose);

%% draw results
if draw
    draw_graphs(simu_desc);
end
%% save information
if saved
    path = pwd;
    path = strcat(path,'\simu_results');
    c = num2cell(clock);
    name = ['Distributed_computing_date_',...
                    num2str(c{1}),'_', num2str(c{2}),'_', num2str(c{3}),'_', num2str(c{4}),'_', num2str(c{5})];
    save([path, '/', name, '.mat'],'simu_desc');
end

end
