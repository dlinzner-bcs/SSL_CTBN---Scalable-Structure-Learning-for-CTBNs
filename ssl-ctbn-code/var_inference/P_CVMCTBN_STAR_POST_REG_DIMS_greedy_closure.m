function [MU,RHO,F,m] = P_CVMCTBN_STAR_POST_REG_DIMS_greedy_closure(node,dt,M,t0,Z,TZ,thresh)
%%%%%%%%%%%%%CLUSTER VARIATIONAL STAR-APPROXIMATION FOR CTBNs WITH ERROR MODEL%%%%%
%%%%%%INPUT:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%node is NODE STRUCT. DEFINING CTBN
%D is local dim. (only D=2 possible a.t.m.)
%T is SIMULATION TIME
%dt is time-step
%M is max. number of iterations for iterative ODE solver
%X0 is initial condition
%Z is noisy observation of state
%TZ is time of observation
%%%%%%OUTPUT:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mu(m,t,i,d) is MARGINAL PROB. IN m'th IERATION OF NODE i TO BE IN STATE d
%AT TIME t
%rho(m,t,i,d) is LAGRANGE MULTIPLIER " " " " " " " " ".
%qt(t,i,d) is AVERAGED RATE """""""
%pst(t,i,d) is AVERAGED CHILD RESPONSE " " " " ".

%%%%%%%%%%%READ OUR SIM. PARAM.%%%%%%%%%%%%%%%%%%

L=length(node);
T=TZ(end)+dt*t0;
tau=ceil(TZ(end)/dt)+t0;
e0=squeeze(Z(:,1,:));


xt=linspace(0,T,tau);
%xt2=linspace(0,T-dt,tau-1);

%%%%%%%%%%SET INITIAL CONDITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%CREATE ANSATZ TRAJECTORY BY PICKING RND CIM%%%
for i=1:L
    D=node(i).D;
    mu_init=zeros(tau,D);
    for d=1:D
        mu_init(1,d)=1/D;
    end
    
    TSPAN = linspace(0,T,tau);
    %%%PICK RANDOM CIM FOR ANSATZ SOLUTION
    q=zeros(tau,D,D);
    Q_i=node(i).cellOfCondRM;
   % c=node(i).allPosStatesOfParents;
   % sc=size(c);
    Q=Q_i{1};
    
    for d=1:D
        for d_=1:D
            q(:,d,d_)=Q(d,d_).*ones(tau,1);
        end
    end
    
    [~,Y] = ode15s(@(t,y) CVM_CTBN_mu_fastD(t, y,D, xt, q(1:tau,:,:), xt, ones(tau,D)), TSPAN,  mu_init(1,:)');
    mu_init(1:tau,:)=Y(1:tau,:);
    for m=1:M
      MU{m}{i}=mu_init(1:tau,:);
      RHO{m}{i}=zeros(tau,D);
    end
end



%%%%%%%%%%ACTUAL SIMULATION OF CTBN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TSPAN_R = linspace(T,0,tau);
TSPAN_M = linspace(0,T,tau);
%options=odeset('RelTol',1e-12,'AbsTol',1e-13);
%%%%%%%%%%%%%%%%%%%%%%%%ITERATE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:M-1
    %%%%%%%%%%UPDATE LAGRANGE MULT. GIVEN MARG. PROBS%%%%%%%%%%%%%%%%%%%%%%
    for i=1:L
        D=node(i).D;
        psi0=zeros(tau,D);
        mu=MU{m}{i};
        rho=RHO{m}{i};
        
        [q_am,q_gm] = P_CVMCTBN_EFF_RATES_SPARSE_REG_DIMS_greedy_approx(node,i,m,MU);
        
        if m>1
             %psi0 = P_CVMCTBN_COUP_D(node,i,squeeze(mu(m,:,:,:)),squeeze(rho(m-1,:,:,:)),D);
            % psi0 = P_CVMCTBN_COUP_DIMS(node,i,m,MU,RHO);
            [psi0] = P_CVMCTBN_COUP_SPARSE_REG_DIMS_greedy(node,i,m,MU,RHO);
        end
      
        rhoT=ones(1,D);
        TSPAN_Y=linspace(T,TZ(end),t0);
        %%%%%%%%%%ACTUAL TIME EVOLUTION OF LAGR. MULT. (BACKWARDS IN TIME)%
        %options = odeset('Jacobian',@(t,y)J_CVM_CTBN_mu(D,t, xt, q(1:tau,:), xt, squeeze(rho(m,1:tau,i,:))));
        [~, Y] = ode15s(@(t,y) CVM_CTBN_rho_fast_sparse_reg_DIMS(D,t, y, xt, q_am, q_am, psi0), TSPAN_Y, rhoT');
        R{length(TZ)+1}=Y;
        
        for k=length(TZ):-1:1
            Y_mem=Y;
            rhoT=Y(end,:).*(squeeze(Z{i}(:,k))')/sum(Y(end,:).*(squeeze(Z{i}(:,k))'));
            if k>1
                TSPAN_Y=linspace(TZ(k),TZ(k-1),ceil((TZ(k)-TZ(k-1))/dt));
            else
                TSPAN_Y=linspace(TZ(1),0,ceil(TZ(1)/dt));
            end
            
            %options = odeset('Jacobian',@(t,y)J_CVM_CTBN_mu(D,t, xt, q(1:tau,:), xt, squeeze(rho(m,1:tau,i,:))));
            [~, Y] = ode15s(@(t,y) CVM_CTBN_rho_fast_sparse_reg_DIMS(D,t, y, xt, q_am,q_gm, psi0), TSPAN_Y, rhoT');
            R{k}=Y;
            [msglast, msgidlast] = lastwarn;
            if isempty(msglast)==0
                Z{i}(:,k)=1; %remove data-point
                warning('') %clear last warning
                msglast=[];
                sprintf('Warning: Could not process data-point %d of node %d',k,i)
                
                rhoT=Y(end,:).*(squeeze(Z{i}(:,k))')/sum(Y(end,:).*(squeeze(Z{i}(:,k))'));
                if k>1
                    TSPAN_Y=linspace(TZ(k),TZ(k-1),ceil((TZ(k)-TZ(k-1))/dt));
                else
                    TSPAN_Y=linspace(TZ(1),0,ceil(TZ(1)/dt));
                end           
                %options = odeset('Jacobian',@(t,y)J_CVM_CTBN_mu(D,t, xt, q(1:tau,:), xt, squeeze(rho(m,1:tau,i,:))));
                [Tt Y] = ode15s(@(t,y) CVM_CTBN_rho_fast_sparse_reg_DIMS(D,t, y, xt, q_am,q_gm, psi0), TSPAN_Y, rhoT');
                R{k}=Y;
            end
        end
        for d=1:D
            A=[];
            for k=length(TZ)+1:-1:1
                A=[A ;R{k}(:,d)];
            end
           % rho(m,tau:-1:1,i,d)=A(1:tau);
            rho(tau:-1:1,d)=A(1:tau);
        end
        RHO{m+1}{i}=rho;
        %%%%%%%%%%UPDATE MARG. PROBS GIVEN LAGRANGE MULT. %%%%%%%%%%%%%%%%%
        %%%%%%%%%%ACTUAL TIME EVOLUTION OF MARG. PROB. (FORWARDS IN TIME)%%
        %  options = odeset('Jacobian',@(t,y)J_CVM_CTBN_mu(D,t, xt, q(1:tau,:), xt, squeeze(rho(m,1:tau,i,:))));
       % [Tt Y] = ode15s(@(t,y) CVM_CTBN_mu_fastD(D,t, y, xt, q(1:tau,:,:), xt, RHO{m}{i}, TSPAN_M, mu(m,1,i,:)));
         [~, Y] = ode15s(@(t,y) CVM_CTBN_mu_fastD(t, y,D, xt, q_gm(1:tau,:,:), xt, rho(1:tau,:)), TSPAN_M, squeeze(mu(1,:)));
     
       %mu(m+1,1:tau,i,:)=Y(1:tau,:);
        MU{m+1}{i}=Y(1:tau,:);
        
        
        if m==M-1
             node(i).R_am=q_am(1:tau,:,:);
        end
        
    end
    
    %%%%%CALCULATE LIKELIHOOD LOWER BOUND TO CHECK FOR CONVERGENCE
%     if m>1
%         [F(m),~,~] = CVM_CTBN_likelihood_starD(node,squeeze(mu(m,:,:,:)),squeeze(rho(m-1,:,:,:)),dt);
%     end
%     if m>5
%         stdF=std(F(end:-1:end-4));
%         if stdF<thresh
%             break
%         end
%         
%     end
    

end
%%%%CONVERGED SOLUTIONS
%MU=squeeze(mu(m,1:tau-t0,:,:));
%RHO=squeeze(rho(m-1,1:tau-t0,:,:));
F=0;
end



