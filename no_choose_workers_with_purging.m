function parameters_out = no_choose_workers_with_purging(parameters)



%for each possible s
s=parameters.s;
t=parameters.t;

%code characteristics
K=t^2*(2*s-1);
Iin = 2*K*(parameters.N)^2/(s*t); %input communication size
Iout = K*(parameters.N)^2/(t^2); %output communication size
I =  Iin+ Iout; %total communication size
C = K*(parameters.N)^3/(s*t^2); %job's number of computation
D_comm_in_vec = Iin./parameters.c_rate_vec_init; %input comm delay
D_comm_out_vec= Iout./parameters.c_rate_vec_init; %output comm delay
Denc =  K*(parameters.N)^2/parameters.mu_enc; %encoding time
Ddec = ((parameters.N)^2* K+ K^3)/parameters.mu_dec; %decoding time
mu_rate_vec_job = parameters.mu_rate_vec_init/C; %service rate of the workers, in job/slot
lambda_J_job = parameters.lambda_J; %task arrival rate
r_comm_scaled =  parameters.c_rate_vec_init/(I*lambda_J_job);%communication rate
r_comp_scaled = mu_rate_vec_job/lambda_J_job;%computation rate
m2_job = 2./mu_rate_vec_job.^2;
a=0.5*lambda_J_job*m2_job.*mu_rate_vec_job;
alpha=parameters.kappa/sum(r_comp_scaled);
phi_l=max(alpha * r_comp_scaled,5/K);%lower bound on utilization

%check hypothesis
OK_index_hyp_phiK = phi_l<=r_comp_scaled;
%1 encoding : mu_enc/N^2 > mu_t
OK_index_enc = 1/Denc >=mu_rate_vec_job;
%2 comm in : c_p/I_in > mu_t
OK_index_comm_in = 1./D_comm_in_vec >=mu_rate_vec_job;
%3 comm out : c_p/I_out > mu_t
OK_index_comm_out = 1./D_comm_out_vec >=mu_rate_vec_job;
%4 decoding : 
OK_index_dec = 1/Ddec  >=mu_rate_vec_job;

Ok_index = logical(OK_index_hyp_phiK.*OK_index_enc.*OK_index_comm_in.*OK_index_comm_out.*OK_index_dec);
r_comp_scaled(~Ok_index)=0;  %computational power of non Ok workers is set to 0


% non_ok=1;
% while non_ok==1
%         %Frow workers satisying the hypothesis, choose subset
    mr_comp_sorted = r_comp_scaled;
    %     
    cum_sum_task = cumsum(mr_comp_sorted);
    Nworkers_tot = find(cum_sum_task>=1+parameters.Theta,1);
    idx_workers_log = mr_comp_sorted(1:Nworkers_tot)>0;
    idx_workers=find(idx_workers_log);
    Nworkers=sum(idx_workers_log);
    
    %compute delay
%     r_comm_scaled_final = r_comm_scaled(idx_workers);
%     r_comp_scaled_final = r_comp_scaled(idx_workers);
%     a_final=a(idx_workers);
%     phi_l_final = phi_l(idx_workers);
%     [phi,~] = optimum_split(r_comp_scaled_final,r_comm_scaled_final,a_final,phi_l_final);
% 
% %     %if some phi do not respect phi_pK>>1, remove workers
%     tol = 5/K;
%     idx_nonzeros = phi>tol;
%     r_comp_scaled(idx_workers(~idx_nonzeros))=0;
%     
%     idx_workers=idx_workers(idx_nonzeros);
%     Nworkers=length(idx_workers);
%     
%     r_comp_scaled_test = mu_rate_vec_job(idx_workers)/lambda_J_job;%computation rate
%     non_ok = sum(r_comp_scaled_test)<1+parameters.Theta;
% end
%select the best s
parameters.mu_rate_vec = parameters.mu_rate_vec_init(idx_workers);
parameters.c_rate_vec = parameters.c_rate_vec_init(idx_workers);
%parameters.var_vec = parameters.var_vec_init(idx_workers);

parameters.P = Nworkers;
%once we know s,t, we can compute the code characteristics
parameters.K=K;
parameters.Iin = Iin;
parameters.Iout = Iout;
parameters.I =  I; %total communication size
parameters.C = C; %task's number of computation
parameters.D_comm_in_vec = D_comm_in_vec(idx_workers); %input comm delay
parameters.D_comm_in_vec_task = parameters.D_comm_in_vec/(parameters.K);
parameters.D_comm_out_vec= D_comm_out_vec(idx_workers); %output comm delay
parameters.D_comm_out_vec_task = parameters.D_comm_out_vec/(parameters.K);
parameters.Denc =  Denc; %encoding time
parameters.Ddec = Ddec; %decoding time
parameters.mu_rate_vec_job = parameters.mu_rate_vec/parameters.C; %service rate of the workers, in task/slot
parameters.mu_rate_vec_task=parameters.mu_rate_vec_job*parameters.K;
parameters.lambda_J_job = parameters.lambda_J; %task arrival rate
parameters.r_comm_scaled =  parameters.c_rate_vec/(parameters.I*parameters.lambda_J_job);%communication rate
parameters.r_comp_scaled = parameters.mu_rate_vec_job/parameters.lambda_J_job;%computation rate
parameters.m2_job = 2./parameters.mu_rate_vec_job.^2;
parameters.a=0.5*parameters.lambda_J_job*parameters.m2_job.*parameters.mu_rate_vec_job;
parameters.alpha=parameters.kappa/sum(parameters.r_comp_scaled);
parameters.phi_l=parameters.alpha * parameters.r_comp_scaled;%lower bound on utilization
 
parameters_out=parameters;
end
