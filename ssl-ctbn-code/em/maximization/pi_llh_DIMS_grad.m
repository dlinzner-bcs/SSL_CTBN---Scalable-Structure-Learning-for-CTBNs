function [F,gF] = pi_llh_DIMS_grad(node,i,lam,pi)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
node(i).pi=pi;
%[node] = ctbn_summarize_stats(node);
%[node] = ctbn_compute_post_rates(node);
F = -marg_llh_sparse_reg_noexpln_single_DIMS(node,i,lam);

if nargout > 1 % gradient required
    D=node(i).D;
    subsets=node(i).subsets;
    for k=1:length(subsets)
        E=0;
        states=node_states(node,subsets{k});
        if isempty(states)==0
            sp=size(states);
            for u=1:sp(1)
                for d=1:D
                    beta_bar=squeeze(node(i).pi(k)*node(i).dwell_k{k}(u,d)+node(i).beta_k{k}(u,d));
                    for d_=1:D
                        if(d~=d_)
                            alpha_bar=squeeze(node(i).pi(k)*node(i).trans_k{k}(u,d,d_)+node(i).alpha_k{k}(u,d,d_));
                            E = E-node(i).trans_k{k}(u,d,d_).*log(beta_bar)-node(i).dwell_k{k}(u,d)*alpha_bar/beta_bar+psi(alpha_bar)*node(i).trans_k{k}(u,d,d_);
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
                        E = E-node(i).trans_k{k}(u,d,d_).*log(beta_bar)-node(i).dwell_k{k}(u,d)*alpha_bar/beta_bar+psi(alpha_bar)*node(i).trans_k{k}(u,d,d_);
                    end
                end
            end
        end
        gF(k) =-( E+lam/node(i).pi(k));
    end
end

end

