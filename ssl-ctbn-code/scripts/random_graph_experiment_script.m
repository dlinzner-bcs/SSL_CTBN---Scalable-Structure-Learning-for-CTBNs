addpath(genpath('./ssl-ctbn-code'))
%Generate synthetic data
L=5; %number of nodes of random graph
max_par=2;
num_graphs=30;
N_TRAJ=10; %number of synthetic trajectories
mworkers=4;
for graphs=1:num_graphs
    A=zeros(L);
    for i=1:L
        lin=[1:L];
        lin(i)=[];
        par=find(mnrnd(1,ones(1,max_par+1).*1/(max_par+1)),1)-1;
        A(i,randsample(lin,par))=1;
    end 
    A
    B=ones(L,L);
    for i=1:L
        B(i,i)=0;
    end
    SIGMA=0.2; %noise in synthetic trajectories
    MAX_PAR=2; %maximum number of parents in synthetic experiments
    steps=10; %number of sampled transitions
    dt=0.005; %simulation time step
    %parameters of ground-truth (kinetic ising-model)
    ta=1; %rate scale
    b=-0.7; % inverse temperature
    D0=2;
    
    states=ones(1,L)*D0;
    %draw samples from ground-truth graph
    node0=createLibOfNodesD(L,states,A, -ones(1,L), b,ta);
    steps0=10; %number of sampled states
    dt0=0.01; %granularity of sample trajectories
    
    %draw sample trajectories
    [time0,DATA0] = CTBN_TRAJ_D(node0,N_TRAJ,steps,steps0,dt0,dt);
    [DATAC,D] = corrupted_observation_gaussianD(DATA0,SIGMA,node0);
    
    name=sprintf('data.mat');
    save(name,'DATAC','time0','L','states');
    
    header=sprintf('states_%d_maxpar_%d_ntraj%d_beta%.2g',D0,max_par,N_TRAJ,b);
    hashA=bi2de(reshape(A,[1,L^2]));
    name=sprintf('ctbn_L%d_hashA%d_b%.2g_tau%.2g',L,hashA,b,ta);
    name=[name '_' header '.mat'];
    %start exhaustive experiment
    ctbn_gradient_structure_learning_dims(name,mworkers)
    %start greedy experiment
    ctbn_gradient_structure_learning_dims_greedy(name,mworkers,2)
end