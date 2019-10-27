function [] = ctbn_gradient_structure_learning_dims_greedy(name,mworkers,K)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load('params.mat');
load('data.mat');

prior.graph=prior_graph;
prior.pi=1;
prior.alpha=alph;
prior.beta=bet;

node=createMotherCTBN_DIMS_K(L,prior,states,K);
[node]=genGammaPriorRatesD_single(node);

delete(gcp('nocreate'))
parpool(mworkers);

for m=1:MAX_SWEEPS
    for k=1:MAX_ITER
        [MU,RHO,node] = ctbn_expectation_sparse_reg_par_DIMS_greedy(node,dt,M,t0,DATAC,time0,thresh);
        f(k)=marg_llh_sparse_reg_DIMS_greedy(node,dt,MU,RHO,DATAC,time0)
    end
    F(m) = marg_llh_sparse_reg_DIMS_greedy(node,dt,MU,RHO,DATAC,time0)
    [node,rep] = estimate_pi_sparse_reg_par_DIMS_grad(node,lam,MAX_RESTARTS);
    
    C=extract_net_prob(node)
    Repm{m}=rep;
    Nodem{m}=node;
    save(name,'node','F','Nodem','Repm')
end


end