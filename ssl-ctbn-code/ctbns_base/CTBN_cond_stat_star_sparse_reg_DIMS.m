function [T,M] = CTBN_cond_stat_star_sparse_reg_DIMS(node,i,t,dt,MU,RHO)
%CALCULATE COND STATISTICS OF CTBN
mu=MU{i};
rho=RHO{i};
D=node(i).D;
%%%get effective rate
ps=node(i).allPosStatesOfParents;
Q_i=node(i).R_gm;
sp=size(ps);
T=zeros(sp(1),D);
M=zeros(sp(1),D,D);
if (isempty(node(i).parents)==0)
    %%%%averaging    
    sp=size(ps);
    for mn=1:sp(1)
        T_i=0;
        M_i=0;       
        pa=ones(t,1);
        ki=1;       
        for k=node(i).parents           
              mu_p=MU{k};
              pa=pa.*mu_p(1:t,ps(mn,ki));  
              ki=ki+1;           
        end      
        for d=1:D            
           T_i(d)=trapz(linspace(0,dt*t,t),mu(1:t,d).*pa);          
           for d_=1:D
               if (d~=d_)
                 M_i(d,d_)=trapz(linspace(0,dt*t,t),pa.*mu(1:t,d).*rho(1:t,d_)./rho(1:t,d)).*Q_i{mn}(d,d_);
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
                M_i(d,d_)=trapz(linspace(0,dt*t,t),mu(1:t,d).*(rho(1:t,d_)./rho(1:t,d))).*Q_i{1}(d,d_);               
             end
        end   
    end
     T(1,:,:)=T_i; 
     M(1,:,:)=M_i; 
end

end

