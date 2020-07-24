function [MU,RHO,node] = ctbn_expectation_sparse_reg_par_DIMS_greedy_closure(node,dt,Mmax,t0,Z,TZ,thresh)
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
        for l=1:length(node(i).pi)
            T{k}{i}{l}=0;
            M{k}{i}{l}=0;
        end
    end
end
parfor k=1:length(TZ)
    try
        [mu,rho,~] = P_CVMCTBN_STAR_POST_REG_DIMS_greedy_closure(node,dt,Mmax,t0,Z{k},TZ{k},thresh);
        MU{k}=mu{Mmax};
        RHO{k}=rho{Mmax-1};
        TK=linspace(0,TZ{k}(end),ceil(TZ{k}(end)/dt));
        for i=1:length(node)
            [T0,M0] = CTBN_cond_stat_star_sparse_reg_DIMS_greedy_approx(node,i,length(TK)-2,dt,MU{k},RHO{k});
            for l=1:length(node(i).pi)
                T{k}{i}{l}=T{k}{i}{l}+T0{l};
                M{k}{i}{l}=M{k}{i}{l}+M0{l};
            end
        end
    catch
        MU{k}=[];
        RHO{k}=[];
        sprintf('Warning: Could not process trajectory %d',k)
    end
end
for i=1:length(node)
    for k=1:length(TZ)
        for l0=1:length(node(i).pi)
            if T{k}{i}{l0}~=0
                st=size(T{k}{i}{l0});
                sm=size(M{k}{i}{l0});
            end
        end
    end
    for l0=1:length(node(i).pi)
        T0=zeros(st(1),st(2));
        M0=zeros(sm(1),sm(2),sm(3));
        for k_=1:length(TZ)
            T0=T0+T{k_}{i}{l0};
            M0=M0+M{k_}{i}{l0};
        end
        node(i).trans_k{l0}=M0;
        node(i).dwell_k{l0}=T0;
    end
end
end

