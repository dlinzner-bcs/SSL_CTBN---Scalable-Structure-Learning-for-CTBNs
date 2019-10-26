function [ C] = extract_net_prob( node )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
L=length(node);
C=zeros(L);
for i=1:L
    for k=1:length(node(i).pi)
        for l=1:L
            if find(node(i).subsets{k}==l)
                C(i,l)= C(i,l)+node(i).pi(k);
            end
        end
    end
end
end

