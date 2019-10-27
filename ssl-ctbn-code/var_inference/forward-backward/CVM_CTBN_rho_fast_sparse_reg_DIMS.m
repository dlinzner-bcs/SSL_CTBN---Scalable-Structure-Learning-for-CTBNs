function dydt = CVM_CTBN_rho_fast_sparse_reg_DIMS(D,t, y, ft, q_am, q_gm, g)

%%% AS ODE IS TIME DEPENDENT ON EFF. RATES WE HAVE TO INTERPOLATE

ind=find(ft>=t,1); %take current time-point fast - but only valid for small timesteps
q_am0=squeeze(q_am(ind,:,:));
q_gm0=squeeze(q_gm(ind,:,:));
g0=squeeze(g(ind,:));
M=q_gm0;
for d=1:D
     M(d,d)=-q_am0(d,d)-g0(d);  
end
dydt = -M*y; % Evalute ODE at times t