function [F] = marg_llh_sparse_reg_DIMS_greedy(node,dt,MU,RHO,Z,TZ)
%node : ctbn
%dt   : time-step of time evolution
%MU   : cell of state marginal probabilities
%   Detailed explanation goes here
F=0;

E=0;
H=0;
for i=1:length(node)
    D=node(i).D;
    subsets=node(i).subsets;
    for k=1:length(subsets)
        if isempty(subsets{k})==0
            states=node_states(node,subsets{k});
            sp=size(states);
            for u=1:sp(1)
                for d=1:D
                    for d_=1:D
                        if(d~=d_)
                            alpha_bar=squeeze(node(i).pi(k)*node(i).trans_k{k}(u,d,d_)+node(i).alpha_k{k}(u,d,d_));
                            beta_bar=squeeze(node(i).pi(k)*node(i).dwell_k{k}(u,d)+node(i).beta_k{k}(u,d));
                            if alpha_bar>170
                                E = E+node(i).alpha_k{k}(u,d,d_)*log(node(i).beta_k{k}(u,d))-stirling_log_gamma(node(i).alpha_k{k}(u,d,d_))-(alpha_bar).*log(beta_bar)+stirling_log_gamma(alpha_bar);
                            else
                                E = E+node(i).alpha_k{k}(u,d,d_)*log(node(i).beta_k{k}(u,d))-log(gamma(node(i).alpha_k{k}(u,d,d_)))-(alpha_bar).*log(beta_bar)+log(gamma(alpha_bar));
                            end
                        end
                    end
                end
            end
        else
            u=1;
            for d=1:D
                for d_=1:D
                    if(d~=d_)
                        alpha_bar=squeeze(node(i).pi(k)*node(i).trans_k{k}(u,d,d_)+node(i).alpha_k{k}(u,d,d_));
                        beta_bar=squeeze(node(i).pi(k)*node(i).dwell_k{k}(u,d)+node(i).beta_k{k}(u,d));
                        if alpha_bar>170
                            E = E+node(i).alpha_k{k}(u,d,d_)*log(node(i).beta_k{k}(u,d))-stirling_log_gamma(node(i).alpha_k{k}(u,d,d_))-(alpha_bar).*log(beta_bar)+stirling_log_gamma(alpha_bar);
                        else
                            E = E+node(i).alpha_k{k}(u,d,d_)*log(node(i).beta_k{k}(u,d))-log(gamma(node(i).alpha_k{k}(u,d,d_)))-(alpha_bar).*log(beta_bar)+log(gamma(alpha_bar));
                        end
                    end
                end
            end
        end
    end
end

F=F+E+H;

end

