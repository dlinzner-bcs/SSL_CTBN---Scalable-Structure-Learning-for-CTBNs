function node=createMotherCTBN_DIMS(L,prior,states)
%input: L - number of vertices (number of nodes)
%       A - adjency matrix (node i is a child of node j if edges(i,j)=1)
%       init - vector of states which consist of 1 and -1
%       temp - temperature
%output: library of structures node{i} with rate matrices encoding glauber
%dynamics
A=prior.graph;
pi_0=prior.pi;
alpha_0=prior.alpha;
beta_0=prior.beta;

%empty lib of nodes
lib(L,'node')=struct('state',0);

for i=1:L
    node(i).parents=[];
end

% for each node create families
for i=1:L
    %vector of parents
    connection=A(i,:);
    firstParent=true;
    for j=1:L
        if connection(j)==1 & firstParent==false
            a=node(i).parents;
            node(i).parents=[a, j];
        end
        if firstParent & connection(j)==1
            node(i).parents=[j];
            firstParent=false;
        end
    end
end
%create children consistent with parents
for i=1:L
    node(i).children=[];
end
for i=1:L
    for j=1:length(node(i).parents)
        node(node(i).parents(j)).children=[node(node(i).parents(j)).children, i];
    end
end

for i=1:L    
    node(i).Omega=[1:states(i)]; 
    node(i).D=length([1:states(i)]);
end
[node] = allPosStatesOfParents(node);

for i=1:L
    D=node(i).D;
    sp=size(node(i).allPosStatesOfParents);
    if sp(1)>0
        node(i).alpha=ones(sp(1),D,D).*alpha_0;
        node(i).beta=ones(sp(1),D).*beta_0;
        node(i).trans=ones(sp(1),D,D).*0;
        node(i).dwell=ones(sp(1),D).*0;
    else
        node(i).alpha=ones(1,D,D).*alpha_0;
        node(i).beta=ones(1,D).*beta_0;
        node(i).trans=ones(1,D,D).*0;
        node(i).dwell=ones(1,D).*0;
    end
end

for i=1:L
    [ P ] = PowerSet(node(i).parents);
    node(i).subsets=P;
    
    sps=size(P);
    node(i).pi=pi_0.*ones(1,sps(2));
    node(i).pi(end)=1;
    node(i).pi=node(i).pi/sum(node(i).pi);
end



end