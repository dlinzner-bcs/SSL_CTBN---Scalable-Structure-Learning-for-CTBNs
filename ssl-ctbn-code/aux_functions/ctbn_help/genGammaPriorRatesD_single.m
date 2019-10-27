function [node]=genGammaPriorRatesD_single(node)
%input: L - number of vertices (number of nodes)
%       edges - adjency matrix (node i is a child of node j if edges(i,j)=1)
%       init - vector of states which consist of 1 and -1 only relevant for sampling
%       tau - prior dwelling-times
%       M - prior transition number
%output: library of structures node{i} with random rate matrices, drawn
%from gamma prior with hyperparams tau and M
% for each node create families
%create conditional rate matrixes
L=length(node);
for w=1:L
    D=node(w).D;
    for d=1:D
        for d_=1:D
            if (d~=d_)
                wD = gamrnd(node(w).alpha(1,d,d_),1/node(w).beta(1,d));
                cellOfCRM{1}(d,d_)=wD;
            end
        end
    end
    for d=1:D
        cellOfCRM{1}(d,d)=-sum(cellOfCRM{1}(d,:));     
    end
    node(w).cellOfCondRM=cellOfCRM;
    clear cellOfCRM;
end

end