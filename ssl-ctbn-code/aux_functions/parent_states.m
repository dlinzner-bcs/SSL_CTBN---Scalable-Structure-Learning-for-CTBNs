function [U] = parent_states(node,i)

%input: node - NonParamCTBN
%       i    - index of initialized node

%output: All possible states of parents A
    
    parentsLength=length(node(n).parents);
    
    %create matrix of all possible states of parents
    U=[1];
    for j=node(n).parents
        
        B=[1:length(node(j).Omega)];
        U=combvec(U,B);
        
    end
    
    
    if parentsLength>=1
        U(1,:)=[];
    else
        U=[];
    end
    
   % node(n).allPosStatesOfParents=U';
    U=U';

end

