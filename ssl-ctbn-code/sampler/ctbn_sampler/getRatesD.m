function [rates, generalRate]=getRatesD(numbOfV,node)
%for each node find its current conditional rate matrix 
%depending on its parent`s states
for w=1:numbOfV
    parentsState = checkParentsState(node, w);
    choosedRateM= chooseConditionalRM(parentsState, node, w);
    node(w).currentCondRM=choosedRateM;
end

%vector and sum of rates depending on state
rates=[];
for w=1:numbOfV
    d=node(w).state;
    rates=[ rates, -node(w).currentCondRM(d,d)];
end
generalRate=sum(rates);
end