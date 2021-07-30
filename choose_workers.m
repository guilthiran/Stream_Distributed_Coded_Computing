function parameters_out = choose_workers(parameters)

Ns = parameters.Ns;

%Initialization
s_vec = linspace(1,parameters.m,Ns);
t_vec=parameters.m./s_vec;
Delay = zeros(1,Ns);
Wkrs = cell(1,Ns);
Nwkrs = zeros(1,Ns);

%for each possible s
for ss = 1:Ns
    s=s_vec(ss);
    t=t_vec(ss);
    
    %code characteristics
    K=t^2*(2*s-1);
    Iin = 2*K*parameters.Omega*(parameters.N)^2/(s*t); %input communication size
    Iout = K*(parameters.N)^2/(t^2); %output communication size
    I =  Iin+ Iout; %total communication size
    C = parameters.Omega*K*(parameters.N)^3/(s*t^2); %task's number of computation
    D_comm_in_vec = Iin./parameters.c_rate_vec_init; %input comm delay
    D_comm_out_vec= Iout./parameters.c_rate_vec_init; %output comm delay
    Denc =  K*parameters.Omega*(parameters.N)^2/parameters.mu_enc; %encoding time
    Ddec = ((parameters.N)^2* K+ K^3)/parameters.mu_dec; %decoding time
    mu_rate_vec_job = parameters.mu_rate_vec_init/C; %service rate of the workers, in task/slot
    lambda_J_job = parameters.lambda_J; %task arrival rate
    r_comm_scaled =  parameters.c_rate_vec_init/(I*lambda_J_job);%communication rate
    r_comp_scaled = mu_rate_vec_job/lambda_J_job;%computation rate
    var_task = parameters.var_vec_init*C^2;
    m2_job = K*var_task+K^2./mu_rate_vec_job.^2;
    a=0.5*lambda_J_job*m2_job.*mu_rate_vec_job;
    alpha=parameters.kappa/sum(r_comp_scaled);
    phi_l=alpha * r_comp_scaled;%lower bound on utilization
    
    %check hypothesis
    %1 encoding : mu_enc/N^2 > mu_t
    OK_index_enc = parameters.mu_enc >=mu_rate_vec_job;
    %2 comm in : c_p/I_in > mu_t
    OK_index_comm_in = 1./D_comm_in_vec >=mu_rate_vec_job;
    %3 comm out : c_p/I_out > mu_t
    OK_index_comm_out = 1./D_comm_out_vec >=mu_rate_vec_job;
    %4 decoding : 
    OK_index_dec = parameters.mu_dec/(parameters.N^2+K^2)  >=mu_rate_vec_job;
    
    Ok_index = logical(OK_index_enc.*OK_index_comm_in.*OK_index_comm_out.*OK_index_dec);
    r_comp_scaled(~Ok_index)=0;  %computational power of non Ok workers is set to 0
    
    %Frow workers satisying the hypothesis, choose subset
    [mr_comp_sorted,Idx] = sort(r_comp_scaled,'descend');
%     mr_comp_sorted=r_comp_scaled;
%     Idx = 1:length(mr_comp_sorted);
%     
    cum_sum_task = cumsum(mr_comp_sorted);
    Nworkers = find(cum_sum_task>=1+parameters.Theta,1);
    if Nworkers>0
        idx_workers = Idx(1:Nworkers);
        r_comm_scaled_final = r_comm_scaled(idx_workers);
        r_comp_scaled_final = r_comp_scaled(idx_workers);
        a_final = a(idx_workers);
        phi_l_final = phi_l(idx_workers);

        %compute delay
        [phi,~] = optimum_split(r_comp_scaled_final,r_comm_scaled_final,a_final,phi_l_final);
        
        %if some phi =0, remove workers
        tol = 1e-6;
        idx_nonzeros = phi>tol;
        idx_workers=idx_workers(idx_nonzeros);
        a_final=a_final(idx_nonzeros);
        r_comp_scaled_final=r_comp_scaled_final(idx_nonzeros);
        r_comm_scaled_final=r_comm_scaled_final(idx_nonzeros);
        
        Delay(ss) = Denc+Ddec...
                    +1/(Nworkers*lambda_J_job)*(sum(a_final'.*phi.^2./(r_comp_scaled_final'-phi))+sum(phi./r_comp_scaled_final'))...
                    +1/Nworkers*sum(phi./(r_comm_scaled_final'*lambda_J_job));
        Wkrs{ss}=idx_workers;
        Nwkrs(ss) = Nworkers;
    else
        Delay(ss) = inf;
        Wkrs{ss}=[];
        Nwkrs(ss) = 0;
    end
end
%select the best s
[~,idx_min]=min(Delay);

%corresponding workers
idx_workers_final = Wkrs{idx_min};%final workers
parameters.mu_rate_vec = parameters.mu_rate_vec_init(idx_workers_final);
parameters.c_rate_vec = parameters.c_rate_vec_init(idx_workers_final);
parameters.var_vec = parameters.var_vec_init(idx_workers_final);

parameters.P = Nwkrs(idx_min);
parameters.s=s_vec(idx_min);
parameters.t=t_vec(idx_min);
%once we know s,t, we can compute the code characteristics
parameters.K=parameters.t^2*(2*parameters.s-1);
parameters.Iin = 2*parameters.K*parameters.Omega*(parameters.N)^2/(parameters.s*parameters.t); %input communication size
parameters.Iout = parameters.K*(parameters.N)^2/(parameters.t^2); %output communication size
parameters.I =  parameters.Iin+ parameters.Iout; %total communication size
parameters.C = parameters.Omega*parameters.K*(parameters.N)^3/(parameters.s*parameters.t^2); %task's number of computation
parameters.D_comm_in_vec = parameters.Iin./parameters.c_rate_vec; %input comm delay
parameters.D_comm_in_vec_task = parameters.D_comm_in_vec/(parameters.K*parameters.Omega);
parameters.D_comm_out_vec= parameters.Iout./parameters.c_rate_vec; %output comm delay
parameters.D_comm_out_vec_task = parameters.D_comm_out_vec/(parameters.K);
parameters.Denc =  parameters.K*parameters.Omega*(parameters.N)^2/parameters.mu_enc; %encoding time
parameters.Ddec = ((parameters.N)^2* parameters.K+ parameters.K^3)/parameters.mu_dec; %decoding time
parameters.mu_rate_vec_job = parameters.mu_rate_vec/parameters.C; %service rate of the workers, in task/slot
parameters.mu_rate_vec_task=parameters.mu_rate_vec_job*parameters.K*parameters.Omega;
parameters.lambda_J_job = parameters.lambda_J; %task arrival rate
parameters.r_comm_scaled =  parameters.c_rate_vec/(parameters.I*parameters.lambda_J_job);%communication rate
parameters.r_comp_scaled = parameters.mu_rate_vec_job/parameters.lambda_J_job;%computation rate
parameters.var_task = parameters.var_vec*parameters.C^2;
parameters.m2_job = parameters.K*parameters.var_task+parameters.K^2./parameters.mu_rate_vec_job.^2;
parameters.a=0.5*parameters.lambda_J_job*parameters.m2_job.*parameters.mu_rate_vec_job;
parameters.alpha=parameters.kappa/sum(parameters.r_comp_scaled);
parameters.phi_l=parameters.alpha * parameters.r_comp_scaled;%lower bound on utilization
 
parameters_out=parameters;
end
