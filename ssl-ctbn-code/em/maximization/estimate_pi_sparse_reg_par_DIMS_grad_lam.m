function [lam] = estimate_pi_sparse_reg_par_DIMS_grad_lam(node,lam,restarts)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(node)
    fun=@(x) pi_llh_DIMS_lam(node,i,x,node(i).pi);
    x0=rand;
    A=[];
    b=[];
    Aeq=[];
    beq=[];
    lb=-10;
    ub=10;
    options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false,'Algorithm','sqp');
    lami = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);
    lam(i) = lami;
end

end

