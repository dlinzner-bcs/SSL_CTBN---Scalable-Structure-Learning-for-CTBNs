%performs exhaustive structure search
%needs a params.mat and a data.mat in same folder (see readme)
%output: Nodem: estimated ctbn at different iterations
%        node : final estimated ctbn at different iterations
%        C    : estimated edge probabilities
%        F    : value of objective function at different iterations
addpath(genpath('../ssl-ctbn-code'))
%name of experiment
name=sprintf('my_test_run_%d_%dparents',1);
%number of workers running in parallel
mworkers=4;

ctbn_gradient_structure_learning_dims(name,mworkers)
