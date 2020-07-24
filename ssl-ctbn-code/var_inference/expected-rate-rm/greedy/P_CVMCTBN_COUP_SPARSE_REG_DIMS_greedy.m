function [psi0] = P_CVMCTBN_COUP_SPARSE_REG_DIMS_greedy(node,i,m,MU,RHO)
%CALCULATE CHILD RESPONSE FOR NODE i
%AT TIME t IN ITERATION m
mu=MU{m}{i};
sm=size(mu);
tau=sm(1);
D=node(i).D;
psi0=zeros(tau,D); % time dimension
%%%%% diamond term
for d=1:D
    %%%%AVERAGING RATES OF CHILDREN%%%%%
    for j=node(i).children
        D_=node(j).D;
        %%%%%REMOVE CONSIDERED PARENT FROM PARENT SET NODE i AS STATE IS FIXED
        n=node(j).parents;
        n(n==i)=[];
        %%%%%LOAD CHILDRENS RATE MATRIXES
        Q_am=zeros(tau,D_,D_);
        %%%%%CALC. EXPECTED RATE MATRIXES
        for d_=1:D_
            for d__=1:D_
                if(d_~=d__)
                    for k=1:length(node(j).subsets)
                        if (isempty(node(j).subsets{k}))==0
                            subsets=node(j).subsets;
                            sub=subsets{k};
                            ind=find(sub==i);
                            % states=de2bi(0:(2^length(sub))-1)+1;
                            states=node_states(node,sub);
                            states(states(:,ind)~=d,:)=[];
                            sub(sub==i)=[];
                            states(:,ind)=[];
                            ss=size(states);
                            for um=1:ss(1)
                                ki=1;
                                pam=ones(tau,1);
                                for h=sub
                                    mu_p=MU{m}{h};
                                    pam=pam.*mu_p(:,states(um,ki));
                                    ki=ki+1;
                                    %  pam=pam.*mu(:,h,states(um,ki));
                                end
                                alpha_bar=squeeze(node(j).pi(k)*node(j).trans_k{k}(um,d_,d__)+node(j).alpha_k{k}(um,d_,d__));
                                beta_bar=squeeze(node(j).pi(k)*node(j).dwell_k{k}(um,d_)+node(j).beta_k{k}(um,d_));
                                Q_am(:,d_,d__)=Q_am(:,d_,d__)+node(j).pi(k).*pam.*alpha_bar/beta_bar;
                            end
                        
                        else
                            Q_am(:,d_,d__)=Q_am(:,d_,d__)+node(j).pi(k)*Q_am(:,d_,d__).*ones(tau,1);
                        end
                    end
                end
            end
        end
        %%%%%CALC. CHILD RESPONSE
        mu_p=MU{m}{j};
        rho_p=RHO{m-1}{j};
        
        C=zeros(tau,1);
        for d_=1:D_
            for d__=1:D_
                if (d_~=d__)
                    a=rho_p(:,d_);
                    a(a==0)=10^(-12);
                    C=C+mu_p(:,d_).*(Q_am(:,d_,d__).*rho_p(:,d__)./a-Q_am(:,d_,d__));                    
                end
            end
        end
      %  psi0(1:tau-1,d)=psi0(1:tau-1,d)+p((j==node(i).children))*C(1:tau-1);%C(tau:-1:1);
      psi0(1:tau-1,d)=psi0(1:tau-1,d)+C(1:tau-1);
    end
end

end

