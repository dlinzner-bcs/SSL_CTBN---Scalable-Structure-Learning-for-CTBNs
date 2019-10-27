# Introduction
This is the companion code to the paper 'Scalable Structure Learning for Continuous-time Bayesian Networks from Incomplete Data' [D.Linzner, M.Schmidt and H Koeppl].

Simulations similar to the ones used in the paper can be run with this Python code.

# Installation
The code has been tested on MacOs 10.14.5.
## Prerequisites
This code has been tested with Matlab 2018b


# Running code
Pre-made scripts for running code can be found in folder "scripts". We provide the files "run_ctbn_ssl.mat" and "run_greedy_ctbn_ssl.mat" for exhaustive and greedy version of our method. Please note that only the greedy version can scale to large number of nodes. Both function expect hyper-parameters and data in a specific format with names "params.mat" and "data.mat". We provide a dummy "params.mat" and "data.mat" for illustration

Further, we provide the random_graph_experiment_script.m that can reconstruct data from Figure 3.

# Hyper-parameters
Hyper-parameters are to be stored in "params.mat" and are loaded automatically during runtime. It needs to contain the following objects:
- M :max number of iterations in Algorithm 1 in inner loop (default 20)
- MAX_ITER: number of iterations in Algorithm 1 in outer loop (default 5)
- MAX_SWEEPS: number of EM iterations (default 3)
- alphabet : shape and rate Gamma params of rate prior
- dt       : time-grid for storing ODE solutions (default dt=0.005)
- lam      : Dirichlet concentration via lam=-(c-1) (default 0.1)
- prior_graph: adjacency matrix of some overcomplete graph (default full graph)
- thesh      : convergence criterium


# Data format
ctbn_ssl method expect likelihoods of latent states and measurement times stores in "data.mat". Data is then loaded automatically during run-time.
The format is:
-time0 : N_TRAJ dimensional cell, where N_TRAJ is number of trajectories. Each cell contains vector of absolute observation time
-DATAC: N_TRAJ dimensional cell, where N_TRAJ is number of trajectories. Each cell contains again a N dimensional cell, where N is the number of nodes. Thus - likelihoods for latent states are stored node wise. These then contain a D dimensional vector, where D is number of states, per observation time - denoting the likelihood for each latent state
         


