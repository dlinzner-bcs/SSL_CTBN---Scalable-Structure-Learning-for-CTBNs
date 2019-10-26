function [node] = ctbn_compute_post_rates_DIMS(node)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(node)
    D=node(i).D;
    clear R
    clear R_gm
    clear R_am
    sp=size(node(i).allPosStatesOfParents);
    if sp(1)>1
        for d=1:sp(1)
            clear R_k
            clear R_k_am
            clear R_k_gm
            R_k=zeros(length(node(i).subsets),D,D);
            state=node(i).allPosStatesOfParents(d,:);
            for k=1:length(node(i).subsets)
                ind=[];
                for l=node(i).subsets{k}
                    ind=[ind,find(node(i).parents==l)];
                end
                states=state(ind);
                j=1;
                for l=1:length(ind)
                    j=j.*(node(i).allPosStatesOfParents(:,ind(l))==states(l))';
                end
                state_ind=find(j==1);
                M_k=squeeze(sum(node(i).trans(state_ind,:,:),1));
                T_k=squeeze(sum(node(i).dwell(state_ind,:),1));
                alpha_k=squeeze(sum(node(i).alpha(state_ind,:,:),1));
                beta_k=squeeze(sum(node(i).beta(state_ind,:),1));
                for d_=1:D
                    for d__=1:D
                        if (d_~=d__)
                            R_k(k,d_,d__)=(node(i).pi(k)*M_k(d_,d__)+ alpha_k(d_,d__))./(node(i).pi(k)*T_k(d_)+beta_k(d_));
                        end
                    end
                end
                R_k_am(k,:,:)=node(i).pi(k).*R_k(k,:,:);
                R_k_gm(k,:,:)=R_k(k,:,:).^node(i).pi(k);
                for h=1:D
                    R_k(k,h,h)=0;
                    R_k_am(k,h,h)=0;
                end
                for h=1:D
                    R_k(k,h,h)=-sum( R_k(k,h,:));
                    R_k_am(k,h,h)=-sum(R_k_am(k,h,:));
                end
                
            end
            R{d}=R_k;
            R_am{d}=squeeze(sum(R_k_am,1));
            R_gm0=squeeze(prod(R_k_gm,1));
            for h=1:D
                R_gm0(h,h)=0;
            end
            for h=1:D
                R_gm0(h,h)=-sum( R_gm0(h,:));
            end
            R_gm{d}= R_gm0;
        end
        
    else
        R{1}=squeeze((node(i).trans(1,:,:)+node(i).alpha(1,:,:))./(node(i).dwell(1,:)+node(i).beta(1,:)));
        R_am{1}=R{1};
        R_gm{1}=R{1};
    end
    node(i).R=R;
    node(i).R_am=R_am;
    node(i).R_gm=R_gm;
end

%evaluate prob j is child of i
for i=1:length(node)
    clear p
    cs=node(i).children;
    if isempty(cs)==0
        for l=1:length(cs)
            prob=0;
            power=node(cs(l)).subsets;
            for k=1:length(power)
                prob=prob+ node(cs(l)).pi(k)*sum(power{k}==i);
            end
            p(l)=prob;
        end
    else
        p=[];
    end
    node(i).child_prob=p;
end


end

