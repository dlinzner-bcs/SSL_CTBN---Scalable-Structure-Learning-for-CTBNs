function [q,logq] = P_CVMCTBN_EFF_RATES_DIMS(node,i,m,MU)
%CALCULATE EFFECTIVE RATES FOR NODE i GIVEN THE NEIGHBOURS
%AT TIME t IN ITERATION m
mu=MU{m}{i};
D=node(i).D;
sm=size(mu);
tau=sm(1);

q=zeros(tau,D,D);
logq=zeros(tau,D,D);

%%%get effective rate
ps=node(i).allPosStatesOfParents;

Q_avg=zeros(tau,D,D);
logQ_avg=zeros(tau,D,D);
Q_i=node(i).cellOfCondRM;
%%%%averaging
for d=1:D
    for d_=1:D
        if (isempty(node(i).parents)==0)
            sp=size(ps);
            for mn=1:sp(1)
                pa=ones(tau,1);
                ki=1;
                
                for k=node(i).parents
                    mu_p=MU{m}{k};
                    pa=pa.*mu_p(:,ps(mn,ki));
                    ki=ki+1;
                end
                
                Q_avg(:,d,d_)=Q_avg(:,d,d_)+Q_i{mn}(d,d_).*pa;
                logQ_avg(:,d,d_)=logQ_avg(:,d,d_)+log(abs(Q_i{mn}(d,d_))).*pa;
            end
            
            %%%%%%%
            
        else
            Q_avg(:,d,d_)=Q_i{1}(d,d_).*ones(tau,1);
            logQ_avg(:,d,d_)=log(abs(Q_i{1}(d,d_)).*ones(tau,1));
        end
    end
end

q=Q_avg;
logq=logQ_avg;

end

