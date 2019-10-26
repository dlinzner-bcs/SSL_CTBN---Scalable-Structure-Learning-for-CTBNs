function [DATAC,D] = corrupted_observation_gaussianD(DATA,SIGMA,node)

N_TRAJ=length(DATA);
sD=size(DATA{1});
L=sD(2);
D=cell(N_TRAJ,1);

for j=1:N_TRAJ
    for i=1:L
        D_h=0;
        dat=squeeze(DATA{j}(:,i));
        D0=node(i).D;
        tau=size(dat);
        Y=zeros(D0,tau(1));
        
        D_h=dat+normrnd(0,ones(tau(1),1)*SIGMA);
        for d_=1:D0
            Y(d_,:)=normpdf(D_h,d_,SIGMA);
        end
        Z=round(Y,5);
        Y(Z==0)=10^(-4);
        Y=Y./sum(Y,1);
        D{j}{i}=D_h;
        DATAC{j}{i}=Y;
    end
end


end

