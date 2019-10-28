%performs greedy structure search
%needs a params.mat and a data.mat in same folder (see readme)
%output: Nodem: estimated ctbn at different iterations
%        node : final estimated ctbn at different iterations
%        C    : estimated edge probabilities
%        F    : value of objective function at different iterations
addpath(genpath('./ssl-ctbn-code'))
%maximal number of parents in greedy search
K=2;
%name of experiment
name=sprintf('my_test_run_%d_%dparents',1,K);
%number of workers running in parallel
mworkers=4;
ctbn_gradient_structure_learning_dims_greedy(name,mworkers,K)
