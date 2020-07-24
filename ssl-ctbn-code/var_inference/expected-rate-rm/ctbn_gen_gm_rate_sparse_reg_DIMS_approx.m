function [q_gm] = ctbn_gen_gm_rate_sparse_reg_DIMS_approx(node,MU,i,m,d,d_)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
mu=MU{m}{i};
sm=size(mu);
tau=sm(1);
subsets=node(i).subsets;
Q_gm=ones(tau,1);
for k=1:length(subsets)
    Q_gm_k=zeros(tau,1);
    sub    = subsets{k};
    states = node_states(node,sub);
    ss     = size(states);
    if isempty(sub)==0
        for um  = 1:ss(1)
            ki  = 1;
            pam = ones(tau,1);
            for j=sub
                mu_p = MU{m}{j};
                pam  = pam.*mu_p(:,states(um,ki));
                ki   = ki+1;
            end
            alpha_bar = node(i).pi(k)*node(i).trans_k{k}(um,d,d_)+node(i).alpha_k{k}(um,d,d_);
            beta_bar  = node(i).pi(k)*node(i).dwell_k{k}(um,d)+node(i).beta_k{k}(um,d);
            Q_gm_k      = Q_gm_k + pam.*(alpha_bar/beta_bar)^node(i).pi(k);
        end
        Q_gm = Q_gm.* Q_gm_k;
    else
        um        = 1;
        pam       = ones(tau,1);
        alpha_bar = node(i).pi(k)*node(i).trans_k{k}(um,d,d_)+node(i).alpha_k{k}(um,d,d_);
        beta_bar  = node(i).pi(k)*node(i).dwell_k{k}(um,d)+node(i).beta_k{k}(um,d);
        Q_gm_k      = Q_gm_k + pam.*(alpha_bar/beta_bar)^node(i).pi(k);
    end
    Q_gm = Q_gm.* Q_gm_k;
end
q_gm = Q_gm;
end

