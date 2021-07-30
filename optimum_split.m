function [phi_middle,overload] = optimum_split(rcomp,rcomm,a,phi_l)
%**************************************************************************
%                       UCLOUVAIN/MIT CONFIDENTIAL
%                           ALL RIGHTS RESERVED
%**************************************************************************

%**************************************************************************
% Function that compute how to split the load between the different
% workers. The arguments are the normalized service rate, the normalized 
% communication rate and the lower bound on phi_p. The output argument gives the
% partition. If the service rate is not sufficient to satisfy the arrival
% rate, the overload variable is set to 1, and the phi are computed by
% increasing artificially the service rate of the worker in order to have a
% stable solution.
%
% Author:           Guillaume Thiran    (guillaume.thiran@uclouvain.be)
%                   2020/12/18
% Last modified :   Guillaume Thiran    (guillaume.thiran@uclouvain.be)
%                   2021/01/26
%**************************************************************************

if ~iscolumn(rcomp)
    rcomp=rcomp';
end
if ~iscolumn(rcomm)
    rcomm=rcomm';
end
if ~iscolumn(a)
    a=a';
end
if ~iscolumn(phi_l)
    phi_l=phi_l';
end
if length(phi_l)==1
    phi_l = phi_l*ones(size(rcomp));
end

overload=0;
%First, check that the service rate is high enough to fulfill the demand
if sum(rcomp)<=1
    %in this case, we choose to artificially increase the service rate to
    %have a solution. Yet, other things could be done
    rcomp=rcomp/sum(rcomp)+eps; 
    overload=1;
end

eta_low = 0; 
fact=10;
eta_high = fact/min(nonzeros(rcomp));
phi_high = compute_phi(rcomp,rcomm,a,phi_l,eta_high);
%check that sum phi_high is greater than 1
iter=0;
iter_max=100;
while sum(phi_high)<=1 && iter<=iter_max
    eta_high=eta_high*2;
    phi_high = compute_phi(rcomp,rcomm,a,phi_l,eta_high);
    iter=iter+1;
end
%At this point we are sure sum(phi_low)<1, sum(phi_high)>1
tol = 1e-9;
iter=0;
iter_max=100;
phi_middle = phi_high;
while abs(sum(phi_middle)-1) >tol && iter < iter_max
    eta_middle = (eta_high+eta_low)/2;
    phi_middle = compute_phi(rcomp,rcomm,a,phi_l,eta_middle);
    
    if sum(phi_middle)>1 %eta is too high
        eta_high = eta_middle;
    else %eta is too low
        eta_low = eta_middle;
    end
    iter=iter+1;
end
end

function phi = compute_phi(rcomp,rcomm,a,phi_l,eta)
    ksi=1./rcomp+1./rcomm-a;
    ok_index = a.*rcomp.^2./((rcomp-phi_l).^2)+ksi <= eta;
    phi_test = rcomp.*(1-sqrt(a./(eta-ksi)));
    phi_test(~ok_index)=phi_l(~ok_index);%if some of them are complex
    phi = max(phi_test,phi_l);
end
