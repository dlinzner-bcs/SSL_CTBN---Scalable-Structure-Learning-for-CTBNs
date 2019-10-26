function [time0,DATA0,DATA,time] = CTBN_TRAJ_D(node,N_TRAJ,steps,steps0,dt0,dt)

L=length(node);

for i=1:N_TRAJ
    tr=0;
    t0=0;
    while (tr<=2*dt)||(t0==0)
        for j=1:L
            ind=randi(length(node(j).Omega));
            node(j).state=node(j).Omega(ind);
        end
        [time0{i},DATA0{i}] = CTBN_Trajectories_IncompleteD(node,steps,steps0,dt0);
        tr=min(diff(time0{i}));
        t0=(min(time0{i}));
        
    end
   
end

end

