function [node,rep] = estimate_pi_sparse_reg_par_DIMS_grad(node,lam,restarts)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
rep=zeros(length(node),restarts);
for i=1:length(node)
    parfor m=1:restarts
        fun=@(x) pi_llh_DIMS_grad(node,i,lam,x);
        x0=rand(1,length(node(i).pi));
        x0=x0/sum(x0);
        sp=length(x0);
        A=[];
        b=[];
        Aeq=ones(1,sp);
        beq=1;
        lb=ones(1,sp)*10^(-3);
        ub=ones(1,sp);  
        options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',true,'Algorithm','sqp');
        argpi = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);
        
        score(m)=fun(argpi);
        PI{m}=argpi;
    end
    [~,ind]=min(score);
    node(i).pi=PI{ind};
    rep(i,:)=score;
    i
end

end

