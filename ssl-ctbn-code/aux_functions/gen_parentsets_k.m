function [S] = gen_parentsets_k(node,i,K)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L=length(node);
x = 1:L;
x(x==i)=[];
for nn = 1:K
    % all the combinations taken nn items at a time
    S0=combnk(x,nn);
    ss=size(S0);
    for i= 1:ss(1)
        if isempty(S0(i,:))==0
          S{nn,i} = S0(i,:);
        end
    end
end
ss=size(S);
S0=reshape(S,[ss(1)*ss(2),1]);
S=[];
k=1;
for i=1:length(S0)
    if isempty(S0{i})==0
    S{k}=S0{i};
    k=k+1;
    end
end
S{k}=[];



end

