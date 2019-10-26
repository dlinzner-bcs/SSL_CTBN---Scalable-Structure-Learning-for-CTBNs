function [] = ctbn_gradient_structure_learning_dims(name,mworkers)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load('params.mat');
load('data.mat');

prior.graph=prior_graph;
prior.pi=1;
prior.alpha=alph;
prior.beta=bet;
node=createMotherCTBN_DIMS(L,prior,states);
for i=1:length(node)
    pi0=ones(length(node(i).pi ),1)*10^(-3);
    pi0(end)=1;
    pi0=pi0/sum(pi0);
    node(i).pi=pi0;
end

[node] = ctbn_summarize_stats_DIMS(node);
[node] = ctbn_compute_post_rates_DIMS(node);

delete(gcp('nocreate'))
parpool(mworkers);

for m=1:MAX_SWEEPS
    %expectation
    for k=1:MAX_ITER
        [MU,RHO,node] = ctbn_expectation_sparse_reg_par_DIMS(node,dt,M,t0,DATAC,time0,thresh);
        f(k)=marg_llh_sparse_reg_DIMS_greedy(node,dt,MU,RHO,DATAC,time0)
    end
    F(m) = marg_llh_sparse_reg_DIMS_greedy(node,dt,MU,RHO,DATAC,time0);
    %maximization
    [node,rep] = estimate_pi_sparse_reg_par_DIMS_grad(node,lam,MAX_RESTARTS);
    C=extract_net_prob(node)
    Repm{m}=rep;
    Nodem{m}=node;
    save(name,'node','F','A','C','Nodem','Repm','pi0','DATA0','time0')
end


end