function node=createLibOfNodesD(numbOfV,states,edges, beginningState, temp,tau)
%input: numbOfV - number of vertices (number of nodes)
%       edges - adjency matrix (node i is a child of node j if edges(i,j)=1)
%       beginning state - vector of states which consist of 1 and -1
%       temp - temperature
%output: libruary of structures node{i}

%empty lib of nodes
lib(numbOfV,'node')=struct('state',0);

for i=1:numbOfV
    node(i).parents=[];
    node(i).Omega=[1:states(i)];
    node(i).D=length(node(i).Omega);
end

% for each node
for i=1:numbOfV
    
    node(i).state=beginningState(i);
    
    %vector of parents
    connection=edges(i,:);
    firstParent=true;
    for j=1:numbOfV
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

[node] = allPosStatesOfParents(node);

%creating conditional rate matrixes
for w=1:numbOfV
    D=node(w).D;
    for d=1:D
        for d_=1:D
            x=((d-d_)-mean(node(w).Omega))/mean(node(w).Omega);
            r=(1/2)*(1+tanh(x*temp*0));
            cellOfCRM{1}(d,d_)=tau*r;
            %cellOfCRM{v}=tau*[-wUp, wUp; wDown, -wDown];
        end
        cellOfCRM{1}(d,d)=0;
        cellOfCRM{1}(d,d)=-sum(cellOfCRM{1}(d,:));
    end
    for v=1:length(node(w).allPosStatesOfParents)
        if length(node(w).parents)>0
            for u=1:length(node(w).parents)
                sumOfParSt=(sum(node(w).allPosStatesOfParents(v,:)-mean(node(w).allPosStatesOfParents,1)));          
            end
        else
            sumOfParSt=0;
        end
        CIM=zeros(D);
        for d=1:D
            for d_=1:D
                if d_~=d
                    x=(d_-mean(node(w).Omega));
                    r=(1/2)*(1+tanh(4*x*temp*sumOfParSt));
                    CIM(d,d_)=tau*r;
                end
            end
            CIM(d,d)=-sum(CIM(d,:));
        end
        cellOfCRM{v}=CIM;
    end
    node(w).cellOfCondRM=cellOfCRM;
    clear cellOfCRM;
end

for i=1:numbOfV
    node(i).children=[];
end

for i=1:numbOfV
    for j=1:length(node(i).parents)
        node(node(i).parents(j)).children=[node(node(i).parents(j)).children, i];
    end
end
end
