function node=createMotherCTBN_DIMS_K(L,prior,states,K)
%input: L - number of vertices (number of nodes)
%       A - adjency matrix (node i is a child of node j if edges(i,j)=1)
%       init - vector of states which consist of 1 and -1
%       temp - temperature
%output: library of structures node{i} with rate matrices encoding glauber
%dynamics
pi_0=prior.pi;
alpha_0=prior.alpha;
beta_0=prior.beta;

%empty lib of nodes
lib(L,'node')=struct('state',0);

for i=1:L
    node(i).parents=[];
end

%create children consistent with parents
for i=1:L
    x=[1:L];
    x(x==i)=[];
    node(i).parents=x;
    node(i).children=x;
end

for i=1:L
    node(i).Omega=[1:states(i)];
    node(i).D=length([1:states(i)]);
end


for i=1:L
    [P] = gen_parentsets_k(node,i,K);
    node(i).subsets=P;
    
    sps=size(P);
    pi=zeros(1,length(P));
    for k=1:length(P)
        if length(P{k})==K
            pi(k)=1;
        end     
    end
    node(i).pi=pi/sum(pi);
end

for i=1:L
    D=node(i).D;
    subsets=node(i).subsets;
    for k=1:length(subsets)
        sp=size(node_states(node,subsets{k}));
        if sp(1)>0
            for u=1:sp(1)
                node(i).alpha_k{k}=ones(u,D,D).*alpha_0/sp(1);
                node(i).beta_k{k}=ones(u,D).*beta_0/sp(1);
                node(i).trans_k{k}=ones(u,D,D).*0;
                node(i).dwell_k{k}=ones(u,D).*0;
            end
        else
            node(i).alpha_k{k}=ones(1,D,D).*alpha_0;
            node(i).beta_k{k}=ones(1,D).*beta_0;
            node(i).trans_k{k}=ones(1,D,D).*0;
            node(i).dwell_k{k}=ones(1,D).*0;
        end
    end
end

for i=1:L
    D=node(i).D;
    node(i).alpha=ones(1,D,D).*alpha_0;
    node(i).beta=ones(1,D).*beta_0;
end


for i=1:length(node)
    clear p
    cs=node(i).children;
    if isempty(cs)==0
        for l=1:length(cs)
            prob=0;
            power=node(cs(l)).subsets;
            for k=1:length(power)
                prob=prob+ node(cs(l)).pi(k)*sum(power{k}==i);
            end
            p(l)=prob;
        end
    else
        p=[];
    end
    node(i).child_prob=p;
end

end