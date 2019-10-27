function [s] = all_node_states(node)
%input: node - NonParamCTBN
%       i    - index of initialized node
%output: All possible states of parents U

U=[1];
for j=1:length(node)
    B=[1:length(node(j).Omega)];
    U=combvec(U,B);
end

U(1,:)=[];

s=U';

end

