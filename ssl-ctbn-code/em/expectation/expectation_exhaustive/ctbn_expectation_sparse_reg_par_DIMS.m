function [MU,RHO,node] = ctbn_expectation_sparse_reg_par_DIMS(node,dt,Mmax,t0,Z,TZ,thresh)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
M=cell(1,length(TZ));
T=cell(1,length(TZ));
MU=cell(1,length(TZ));
RHO=cell(1,length(TZ));

for k=1:length(TZ)
    MU{k}=[];
    RHO{k}=[];
    for i=1:length(node)
        T{k}{i}=0;
        M{k}{i}=0;
    end
end
parfor k=1:length(TZ)
    try
        [mu,rho,~] = P_CVMCTBN_STAR_POST_REG_DIMS(node,dt,Mmax,t0,Z{k},TZ{k},thresh);
        MU{k}=mu{Mmax};
        RHO{k}=rho{Mmax-1};
        TK=linspace(0,TZ{k}(end),ceil(TZ{k}(end)/dt));
        for i=1:length(node)
            [T0,M0] = CTBN_cond_stat_star_sparse_reg_DIMS(node,i,length(TK)-2,dt,MU{k},RHO{k});
            T{k}{i}=T0;
            M{k}{i}=M0;
        end
    catch
        MU{k}=[];
        RHO{k}=[];
        sprintf('Warning: Could not process trajectory %d',k)
    end
end
for i=1:length(node)
    for k=1:length(TZ)
        if T{k}{i}~=0
          st=size(T{k}{i});
          sm=size(M{k}{i});
        end
    end
    T0=zeros(st(1),st(2));
    M0=zeros(sm(1),sm(2),sm(3));
    for k=1:length(TZ)
        T0=T0+T{k}{i};
        M0=M0+M{k}{i};
    end
    node(i).trans=M0;
    node(i).dwell=T0;
end
[node] = ctbn_summarize_stats_DIMS(node);
[node] = ctbn_compute_post_rates_DIMS(node);
end

