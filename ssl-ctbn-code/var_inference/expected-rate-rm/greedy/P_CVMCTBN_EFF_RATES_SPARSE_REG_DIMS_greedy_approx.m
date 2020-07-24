function [q_am,q_gm] = P_CVMCTBN_EFF_RATES_SPARSE_REG_DIMS_greedy_approx(node,i,m,MU)
%CALCULATE EFFECTIVE RATES FOR NODE i GIVEN THE NEIGHBOURS
%AT TIME t IN ITERATION m
mu=MU{m}{i};
D=node(i).D;
sm=size(mu);
tau=sm(1);

%%%get effective rate
Q_am=zeros(tau,D,D);
Q_gm=zeros(tau,D,D);
%%%%averaging
for d=1:D
    for d_=1:D
        if(d~=d_)
            Q_am(:,d,d_) = ctbn_gen_am_rate_sparse_reg_DIMS(node,MU,i,m,d,d_);
            Q_gm(:,d,d_) = ctbn_gen_gm_rate_sparse_reg_DIMS_approx(node,MU,i,m,d,d_);
        end
    end
end
q_am=Q_am;
q_gm=Q_gm;



