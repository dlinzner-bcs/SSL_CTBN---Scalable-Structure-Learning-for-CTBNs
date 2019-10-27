function [q_am] = ctbn_gen_am_rate_sparse_reg_DIMS(node,MU,i,m,d,d_)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
mu=MU{m}{i};
sm=size(mu);
tau=sm(1);
subsets=node(i).subsets;
Q_am=zeros(tau,1);
for k=1:length(subsets)
    sub=subsets{k};
    %states=de2bi(0:(2^length(sub))-1)+1;
    states=node_states(node,sub);
    ss=size(states);
    if isempty(sub)==0
    for um=1:ss(1)
        ki=1;
        pam=ones(tau,1);
        for j=sub
             mu_p=MU{m}{j};
             pam=pam.*mu_p(:,states(um,ki));
           % pam=pam.*mu(:,j,states(um,ki));
             ki=ki+1;
        end
        alpha_bar=squeeze(node(i).pi(k)*node(i).trans_k{k}(um,d,d_)+node(i).alpha_k{k}(um,d,d_));
        beta_bar=squeeze(node(i).pi(k)*node(i).dwell_k{k}(um,d)+node(i).beta_k{k}(um,d));
        Q_am=Q_am+node(i).pi(k).*pam.*alpha_bar/beta_bar;
    end
    else
        um=1;
        pam=ones(tau,1);
        alpha_bar=squeeze(node(i).pi(k)*node(i).trans_k{k}(um,d,d_)+node(i).alpha_k{k}(um,d,d_));
        beta_bar=squeeze(node(i).pi(k)*node(i).dwell_k{k}(um,d)+node(i).beta_k{k}(um,d));
        Q_am=Q_am+node(i).pi(k).*pam.*alpha_bar/beta_bar;
    end
end
q_am=Q_am;
end

