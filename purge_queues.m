function [comm_in,Queue_comm_in,comp,Queue_comp,comm_out,Queue_comm_out] = purge_queues(comm_in,Queue_comm_in,comp,Queue_comp,comm_out,Queue_comm_out,j,P,task_desc)

for p=1:P
    %comm in 
    m = comm_in(p,1); 
    if m>=0 && task_desc(m).job ==j
        comm_in(p,:)=-1;
    end
        
    %Queue_comm_in
    Lin = length(Queue_comm_in{p});
    idx=1;
    while idx <=Lin 
        m=Queue_comm_in{p}(idx);
        if task_desc(m).job ==j
            Queue_comm_in{p}(idx)=[];
            Lin=Lin-1;
        else
            idx=idx+1;
        end
    end
    
    %cop 
    m = comp(p,1); 
    if m>=0 && task_desc(m).job ==j
        comp(p,:)=-1;
    end
        
    %Queue_comp
    Lcomp = length(Queue_comp{p});
    idx =1;
    while idx <=Lcomp 
        m=Queue_comp{p}(idx);
        if task_desc(m).job ==j
            Queue_comp{p}(idx)=[];
            Lcomp=Lcomp-1;
        else
            idx=idx+1;
        end
    end
    
    %comm out 
%     m = comm_out(p,1); 
%     if m>=0 && task_desc(m).job ==j
%         comm_out(p,:)=-1;
%     end
%         
%     %Queue_comm_out
%     Lout = length(Queue_comm_out{p});
%     idx = 1;
%     while idx <=Lout 
%         m=Queue_comm_out{p}(idx);
%         if task_desc(m).job ==j
%             Queue_comm_out{p}(idx)=[];
%             Lout=Lout-1;
%         else
%             idx=idx+1;
%         end
%     end
    
end



end