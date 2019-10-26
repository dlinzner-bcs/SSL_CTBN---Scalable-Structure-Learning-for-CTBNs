function [U] = node_states(node,ind)
%input: node - NonParamCTBN
%       i    - index of initialized node
%output: All possible states of parents U
    
    parentsLength=length(ind);
    %create matrix of all possible states of parents
    U=[1];
    for j=ind        
        B=[1:length(node(j).Omega)];
        U=combvec(U,B);       
    end
    if parentsLength>=1
        U(1,:)=[];
    else
        U=[];
    end
    U=U';

end

