function [F] = marg_llh_sparse_reg_noexpln_single_DIMS(node,i,lam)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
F=0;
E=0;

D=node(i).D;

subsets=node(i).subsets;
for k=1:length(subsets)
    states=node_states(node,subsets{k});
    if isempty(states)==0
        sp=size(states);
        for u=1:sp(1)
            for d=1:D
                beta_bar=squeeze(node(i).pi(k)*node(i).dwell_k{k}(u,d)+node(i).beta_k{k}(u,d));
                for d_=1:D
                    if(d~=d_)
                        alpha_bar=squeeze(node(i).pi(k)*node(i).trans_k{k}(u,d,d_)+node(i).alpha_k{k}(u,d,d_));
                         E = E-(alpha_bar).*log(beta_bar)+gammaln(alpha_bar);
                    end
                end
            end
        end
    else
        u=1;
        for d=1:D
            beta_bar=squeeze(node(i).pi(k)*node(i).dwell_k{k}(u,d)+node(i).beta_k{k}(u,d));
            for d_=1:D
                if(d~=d_)
                    alpha_bar=squeeze(node(i).pi(k)*node(i).trans_k{k}(u,d,d_)+node(i).alpha_k{k}(u,d,d_));
                    
                        E = E-(alpha_bar).*log(beta_bar)+gammaln(alpha_bar);
                    
                end
            end
        end
    end
end
P=lam*sum(log(node(i).pi));

F=F+E+P;

end

