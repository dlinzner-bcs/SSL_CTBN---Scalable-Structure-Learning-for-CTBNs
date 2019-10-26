function [T_k,M_k] = CTBN_cond_stat_star_sparse_reg_DIMS_greedy(node,i,t,dt,MU,RHO)
%CALCULATE COND STATISTICS OF CTBN
mu=MU{i};
rho=RHO{i};
D=node(i).D;
%%%get effective rate
for k=1:length(node(i).pi)
    ps=node_states(node,node(i).subsets{k});
    sp=size(ps);
    T=zeros(sp(1),D);
    M=zeros(sp(1),D,D);
    if (isempty(node(i).subsets{k})==0)
        %%%%averaging
        sp=size(ps);
        for mn=1:sp(1)
            
            T_i=0;
            M_i=0;
            pa=ones(t,1);
            ki=1;
            for l=node(i).subsets{k}
                mu_p=MU{l};
                pa=pa.*mu_p(1:t,ps(mn,ki));              
                ki=ki+1;
            end
            for d=1:D
                T_i(d)=trapz(linspace(0,dt*t,t),mu(1:t,d).*pa);
                for d_=1:D
                    if (d~=d_)
                        q_am=ctbn_gen_am_rate_sparse_reg_DIMS_partial(node,MU,i,ps(mn,:),node(i).subsets{k},d,d_);
                        M_i(d,d_)=trapz(linspace(0,dt*t,t),q_am(1:t).*pa(1:t).*mu(1:t,d).*(rho(1:t,d_))./(rho(1:t,d)));
                    end
                end
            end
            T(mn,:)=T_i(:);
            M(mn,:,:)=M_i(:,:);
        end
    else
        for d=1:D
            T_i(d)=trapz(linspace(0,dt*t,t),mu(1:t,d));
            for d_=1:D
                if (d~=d_)
                    q_am=ctbn_gen_am_rate_sparse_reg_DIMS_partial(node,MU,i,[],[],d,d_);
                    M_i(d,d_)=trapz(linspace(0,dt*t,t),mu(1:t,d).*q_am(1:t).*(rho(1:t,d_)./rho(1:t,d)));
                end
            end
        end
        T(1,:,:)=T_i;
        M(1,:,:)=M_i;
    end
    T_k{k}=T;
    M_k{k}=M;
end

