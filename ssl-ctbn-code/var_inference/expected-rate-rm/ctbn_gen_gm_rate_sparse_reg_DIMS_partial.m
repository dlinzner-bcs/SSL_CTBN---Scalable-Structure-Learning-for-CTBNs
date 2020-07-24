function [q_gm] = ctbn_gen_gm_rate_sparse_reg_DIMS_partial(node,MU,i,state,par,d,d_)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
mu=MU{i};
sm=size(mu);
tau=sm(1);
subsets=node(i).subsets;
Q_gm=ones(tau,1);
for k=1:length(subsets)
    Q_gm_k=zeros(tau,1);
    sub=subsets{k};
    sub0=sub;
    ind=[];
    ind_l=[];
    for l=1:length(par)
        ind=[ind find(sub==par(l))];
        if isempty(find(sub==par(l)))==0
            ind_l=[ind_l,l];
        end
    end
    
    if isempty(ind)
        ind=[];
        ind_l=[];
    end
    sub0(ind)=[];
    states=node_states(node,sub);
    ss=size(states);
    states0=states;
    states0(:,ind)=[];
    ind_u=ones(ss(1),1);
    for um=1:ss(1)
        if isempty(ind)==0
            ind_u(um)=prod(states(um,ind)==state(ind_l));
        end
    end
    if ss(1)>0
        for um=1:ss(1)
            ki=1;
            pam=ones(tau,1);
            for j=sub0
                mu_p=MU{j};
                pam=pam.*mu_p(:,states0(um,ki));
                ki=ki+1;
            end
            alpha_bar=(squeeze(node(i).pi(k)*node(i).trans_k{k}(um,d,d_)+node(i).alpha_k{k}(um,d,d_)));
            beta_bar=(squeeze(node(i).pi(k)*node(i).dwell_k{k}(um,d)+node(i).beta_k{k}(um,d)));           
            Q_gm_k=Q_gm_k+ind_u(um).*pam.*(alpha_bar/beta_bar)^node(i).pi(k);            
        end
    else
        um=1;
        pam=ones(tau,1);
        alpha_bar=squeeze(node(i).pi(k)*node(i).trans_k{k}(um,d,d_)+node(i).alpha_k{k}(um,d,d_));
        beta_bar=squeeze(node(i).pi(k)*node(i).dwell_k{k}(um,d)+node(i).beta_k{k}(um,d));       
        Q_gm_k=Q_gm_k+pam.*(alpha_bar/beta_bar)^node(i).pi(k);  
    end
    Q_gm=Q_gm_k.*Q_gm;
end
q_gm=Q_gm;
end

