clearvars;
close all;

PP = 150;
Ns = 100;

%% ------------------ Parameters -------------------------------- %
parameters.lambda_J=0.001; %arrival rate of jobs, in [job/slots]
parameters.Omega = 1.5; %redundancy factor
parameters.N = 100; %size of the matrix
parameters.m=50; %number of partitions
parameters.mu_rate_vec_init=1000*rand(1,PP); %service rate of the workers, in [comp/slots]
parameters.c_rate_vec_init = 200*rand(1,PP); %communication rate of the workers, in [bits/slots]
parameters.mu_enc = 10000; %encoding service rate in [comp/slots]
parameters.mu_dec = 100000; %decoding service rate in [comp/slots]
parameters.kappa = 0.2; %lower bound on utilization, sum phi_l = kappa
parameters.Theta =2;

%------------------------Comp/Comme trade-off------------------------%
s_vec = linspace(1,parameters.m,Ns);
t_vec=parameters.m./s_vec;
Delay = zeros(1,Ns);
Delay_enc = zeros(1,Ns);
Delay_comm_in = zeros(1,Ns);
Delay_comp = zeros(1,Ns);
Delay_comm_out = zeros(1,Ns);
Delay_dec = zeros(1,Ns);
Fail_enc= zeros(1,Ns);
Fail_comm_in = zeros(1,Ns);
Fail_comm_out = zeros(1,Ns);
Fail_dec = zeros(1,Ns);
Fail_tot = zeros(1,Ns);

Sum_rcomscaled = zeros(1,Ns);

Wkrs = cell(1,Ns);
Nwkrs = zeros(1,Ns);
for ss = 1:Ns
    s=s_vec(ss);
    t=t_vec(ss);
    
    %once we know s,t, we can compute the code characteristics
    K=t^2*(2*s-1);
    Iin = 2*K*parameters.Omega*(parameters.N)^2/(s*t); %input communication size
    Iout = K*parameters.Omega*(parameters.N)^2/(t^2); %output communication size
    I =  Iin+ Iout; %total communication size
    C = parameters.Omega*K*(parameters.N)^3/(s*t^2); %job's number of computation
    D_comm_in_vec = Iin./parameters.c_rate_vec_init; %input comm delay
    D_comm_out_vec= Iout./parameters.c_rate_vec_init; %output comm delay
    Denc =  K*parameters.Omega*(parameters.N)^2/parameters.mu_enc; %encoding time
    Ddec = ((parameters.N)^2* K+ K^3)/parameters.mu_dec; %decoding time
    mu_rate_vec_job = parameters.mu_rate_vec_init/C; %service rate of the workers, in job/slot
    lambda_J_job = parameters.lambda_J; %task arrival rate
    r_comm_scaled =  parameters.c_rate_vec_init/(I*lambda_J_job);%communication rate
    r_comp_scaled = mu_rate_vec_job/lambda_J_job;%computation rate
    m2_job = 2./mu_rate_vec_job.^2;
    a=0.5*lambda_J_job*m2_job.*mu_rate_vec_job;
    alpha=parameters.kappa/sum(r_comp_scaled);
    phi_l=alpha * r_comp_scaled;%lower bound on utilization
    
    %check hypothesis
    %1 encoding : mu_enc/N^2 > mu_t
    OK_index_enc = 1/Denc >=mu_rate_vec_job;
    %2 comm in : c_p/I_in > mu_t
    OK_index_comm_in = 1./D_comm_in_vec >=mu_rate_vec_job;
    %3 comm out : c_p/I_out > mu_t
    OK_index_comm_out = 1./D_comm_out_vec >=mu_rate_vec_job;
    %4 decoding : 
    OK_index_dec = 1/Ddec  >=mu_rate_vec_job;
    
    Fail_enc(ss)= sum(OK_index_enc);
    Fail_comm_in(ss) = sum(OK_index_comm_in);
    Fail_comm_out(ss) = sum(OK_index_comm_out);
    Fail_dec(ss) = sum(OK_index_dec);
    
    
    Ok_index = logical(OK_index_enc.*OK_index_comm_in.*OK_index_comm_out.*OK_index_dec);
    Fail_tot(ss)=sum(Ok_index);
    r_comp_scaled(~Ok_index)=0;  %computational power of non Ok workers is set to 0
    
    %Frow workers satisying the hypothesis, choose subset
    [mr_comp_sorted,Idx] = sort(r_comp_scaled,'descend');
    cum_sum_task = cumsum(mr_comp_sorted);
    Nworkers = find(cum_sum_task>=1+parameters.Theta,1);
    if Nworkers>0
        idx_workers = Idx(1:Nworkers);
        %mu_final = mu_rate_vec_task(idx_workers);
        c_final = parameters.c_rate_vec_init(idx_workers);
        r_comm_scaled_final = r_comm_scaled(idx_workers);
        r_comp_scaled_final = r_comp_scaled(idx_workers);
        phi_l_final = phi_l(idx_workers);
        a_final=a(idx_workers);

        %compute delay
        [phi,~] = optimum_split(r_comp_scaled_final,r_comm_scaled_final,a_final,phi_l_final);
        Delay_enc(ss)=Denc;
        Delay_comm_in(ss)=1/Nworkers*sum(phi*Iin./c_final');
        Delay_comm_out(ss)=1/Nworkers*sum(phi*Iout./c_final');
        Delay_comp(ss)=1/(Nworkers*lambda_J_job)*sum(a_final'.*phi.^2./(r_comp_scaled_final'-phi)+phi./r_comp_scaled_final');
        Delay_dec(ss)=Ddec;
        Delay(ss) = Denc+Ddec...
                    +Delay_comm_in(ss) + Delay_comm_out(ss)...
                    +Delay_comp(ss);
        Sum_rcomscaled(ss)=sum(r_comp_scaled_final);
        Wkrs{ss}=idx_workers;
        Nwkrs(ss) = Nworkers;
    else
        %Fail_comp(ss)=1;
        Delay_enc(ss)=inf;
        Delay_comm_in(ss)=inf;
        Delay_comp(ss)=inf;
        Delay_comm_out(ss)=inf;
        Delay_dec(ss)=inf;
        Delay(ss) = inf;
        Wkrs{ss}=[];
        Nwkrs(ss) = 0;
        Sum_rcomscaled(ss)=inf;
    end

end

%% graphs
gamma = 0; 
Delay= Delay+gamma*Nwkrs;
figure;
subplot(132);
plot(s_vec,Delay,'b','linewidth',1.5);
grid on;
[Dmin,idx_min]=min(Delay);
hold on;
scatter(s_vec(idx_min),Dmin,50,'r');
title('Delay')
plot(s_vec,Delay_enc,'r','linewidth',1.5);
plot(s_vec,Delay_comm_in,'g','linewidth',1.5);
plot(s_vec,Delay_comp,'k','linewidth',1.5);
plot(s_vec,Delay_comm_out,'m','linewidth',1.5);
plot(s_vec,Delay_dec,'c','linewidth',1.5);
legend('Tot','optimal s','Enc','Comm in', 'Comp', 'Comm out', 'Dec');
xlabel('s')
xlim([0 parameters.m])

subplot(131)
plot(s_vec, Fail_enc,'r','linewidth',1.5);
hold on;
plot(s_vec,Fail_comm_in,'g','linewidth',1.5);
%plot(s_vec,Fail_comp,'k','linewidth',1.5);
plot(s_vec,Fail_comm_out,'m','linewidth',1.5);
plot(s_vec,Fail_dec,'c','linewidth',1.5);
plot(s_vec,Fail_tot,'--b','linewidth',1.5);
title('Number of Ok workers - hypothesis');
xlabel('s');
legend('Enc','Comm in', 'Comm out', 'dec', 'tot');

subplot(133);
plot(s_vec,Sum_rcomscaled,'b','linewidth',1.5);
grid on;
hold on;
scatter(s_vec(idx_min),Sum_rcomscaled(idx_min),50,'r');
title('\Sigma_p r_p^{comp}')
ylim([0 1+parameters.Theta+0.3])
xlim([0 parameters.m])
xlabel('s')

