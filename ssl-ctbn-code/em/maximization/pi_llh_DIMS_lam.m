function [F] = pi_llh_DIMS_lam(node,i,lam,pi)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
node(i).pi=pi;
%[node] = ctbn_summarize_stats(node);
%[node] = ctbn_compute_post_rates(node);
F = -marg_llh_sparse_reg_noexpln_single_DIMS(node,i,lam);

end

