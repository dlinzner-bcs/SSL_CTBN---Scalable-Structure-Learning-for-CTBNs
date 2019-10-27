function [q_am,q_gm] = P_CVMCTBN_EFF_RATES_SPARSE_REG_DIMS(node,i,m,MU)
%CALCULATE EFFECTIVE RATES FOR NODE i GIVEN THE NEIGHBOURS
%AT TIME t IN ITERATION m
mu=MU{m}{i};
D=node(i).D;
sm=size(mu);
tau=sm(1);

%%%get effective rate
ps=node(i).allPosStatesOfParents;

Q_am=zeros(tau,D,D);
Q_gm=zeros(tau,D,D);
Q_i_am=node(i).R_am;
Q_i_gm=node(i).R_gm;
%%%%averaging
for d=1:D
    for d_=1:D
        if(d~=d_)
            if (isempty(node(i).parents)==0)
                sp=size(ps);
                for mn=1:sp(1)
                    pa=ones(tau,1);
                    ki=1;
                    for k=node(i).parents
                        mu_p=MU{m}{k};
                        pa=pa.*mu_p(:,ps(mn,ki));
                        % pa=pa.*mu(:,k,ps(mn,ki));
                        ki=ki+1;
                    end
                    Q_gm(:,d,d_)=Q_gm(:,d,d_)+Q_i_gm{mn}(d,d_).*pa;
                end
                Q_am(:,d,d_) = ctbn_gen_am_rate_sparse_reg_DIMS(node,MU,i,m,d,d_);
            else
                Q_am(:,d,d_)=Q_i_am{1}(d,d_).*ones(tau,1);
                Q_gm(:,d,d_)=Q_i_gm{1}(d,d_).*ones(tau,1);
            end
        end
    end
end
q_am=Q_am;
q_gm=Q_gm;
end

