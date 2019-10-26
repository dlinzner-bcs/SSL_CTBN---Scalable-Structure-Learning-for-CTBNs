function [node] = ctbn_summarize_stats_DIMS(node)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(node)
    D=node(i).D;
    clear alpha
    clear beta
    clear M
    clear T
    if length(node(i).subsets)>1    
        for k=1:length(node(i).subsets)
            sp=size(node(i).allPosStatesOfParents);
            if k==1
                M_k=zeros(1,D,D);
                alpha_k=M_k;
                T_k=zeros(1,D);
                beta_k=T_k;
            
                M_k(1,:,:)=sum(node(i).trans,1);
                T_k(1,:)=sum(node(i).dwell,1);
                alpha_k(1,:,:)=sum(node(i).alpha,1);
                beta_k(1,:)=sum(node(i).beta,1);
            else
                ind=[];
                for l=node(i).subsets{k}
                    ind=[ind,find(node(i).parents==l)];
                end
                %states=de2bi(0:(2^length(ind))-1)+1;
                states=node_states(node,node(i).subsets{k});
                ss=size(states);
                M_k=zeros(ss(1),D,D);
                alpha_k=M_k;
                T_k=zeros(ss(1),D);
                beta_k=T_k;
                for u=1:ss(1)
                    j=1;
                    for l=1:length(ind)
                        j=j.*(node(i).allPosStatesOfParents(:,ind(l))==states(u,l))';
                    end
                    state_ind=find(j==1);                   
                     M_k(u,:,:)=sum(node(i).trans(state_ind,:,:),1);
                     T_k(u,:)=sum(node(i).dwell(state_ind,:),1);
                     alpha_k(u,:,:)=sum(node(i).alpha(state_ind,:,:),1);
                     beta_k(u,:)=sum(node(i).beta(state_ind,:),1);        
                end
            end
            M{k}=M_k;
            T{k}=T_k;
            alpha{k}=alpha_k;
            beta{k}=beta_k;
        end
        
    else
        M_k=zeros(1,D,D);
        alpha_k=M_k;
        T_k=zeros(1,D);
        beta_k=T_k;
        state_ind=1;
        M_k(1,:,:)=sum(node(i).trans(state_ind,:,:),1);
        T_k(1,:)=sum(node(i).dwell(state_ind,:),1);
        alpha_k(1,:,:)=sum(node(i).alpha(state_ind,:,:),1);
        beta_k(1,:)=sum(node(i).beta(state_ind,:),1);
        M{1}=M_k;
        T{1}=T_k;
        alpha{1}=alpha_k;
        beta{1}=beta_k;
    end
    node(i).trans_k=M;
    node(i).dwell_k=T;
    node(i).alpha_k=alpha;
    node(i).beta_k=beta;
end

end

