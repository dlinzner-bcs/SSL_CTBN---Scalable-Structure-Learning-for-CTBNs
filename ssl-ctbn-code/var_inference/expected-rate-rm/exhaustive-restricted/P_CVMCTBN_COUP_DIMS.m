function [psi0] = P_CVMCTBN_COUP_DIMS(node,i,m,MU,RHO)
%CALCULATE CHILD RESPONSE FOR NODE i
%AT TIME t IN ITERATION m
mu=MU{m}{i};
sm=size(mu);
tau=sm(1);
D=node(i).D;
psi0=zeros(tau,D,D);
%%%%% diamond term
for d=1:D
    %%%%AVERAGING RATES OF CHILDREN%%%%%
    for j=node(i).children
       
        D_=node(j).D;
        %%%%%REMOVE CONSIDERED PARENT FROM PARENT SET NODE i AS STATE IS FIXED
        n=node(j).parents;
        n(n==i)=[];
        %%%%%STATE TRANSFORM ALWAYS NECESSARY
        ps=node(j).allPosStatesOfParents;
        %%%%%CONDITION ON STATE OF NODE i BY REMOVING ALL NON CONSITENT STATES
        sp=size(ps);
        pk=[1:sp(1)];
        pk(ps(:,i==node(j).parents)~=d)=[];
        ps(ps(:,i==node(j).parents)~=d,:)=[];
        ps(:,i==node(j).parents)=[];
        %%%%%LOAD CHILDRENS RATE MATRIXES
        Q_i=node(j).cellOfCondRM;
        Q_avg=zeros(tau,D_,D_);
        sp=size(ps);
        %%%%%CALC. EXPECTED RATE MATRIXES
        for d_=1:D_
            for d__=1:D_
                if (isempty(node(j).parents))==0
                    for mn=1:sp(1)
                        pa=ones(tau,1);
                        ki=1;
                        for k=n
                            mu_p=MU{m}{k};
                            pa=pa.*mu_p(:,ps(mn,ki));
                            ki=ki+1;
                        end
                        Q_avg(:,d_,d__)=Q_avg(:,d_,d__)+ Q_i{pk(mn)}(d_,d__).*pa;
                    end
                else
                    Q_avg(:,d_,d__)=Q_i{1}(d_,d__).*ones(tau,1);
                end
            end
        end
        %%%%%CALC. CHILD RESPONSE
        C=zeros(tau,1);
        mu_p=MU{m}{j};
        rho_p=RHO{m}{j};
        for d_=1:D_
            if prod(rho_p(:,d_)>0)
                for d__=1:D_
                    if d_~=d__
                        %  if prod(rho_p(:,d_)>0)
                        C=C+mu_p(:,d_).*(Q_avg(:,d_,d__).*rho_p(:,d__)./rho_p(:,d_));
                        %   end
                    end
                end
                C=C+mu_p(:,d_).*Q_avg(:,d_,d_);
            end
        end
        psi0(1:tau-1,d)=psi0(1:tau-1,d)+C(1:tau-1);%C(tau:-1:1);
    end
end
end



