function dydt = CVM_CTBN_rho_fastD(D,t, y, ft, f, gt, g)

%%% AS ODE IS TIME DEPENDENT ON EFF. RATES WE HAVE TO INTERPOLATE


% for d=1:D
%     for d_=1:D
%         f0(d,d_)=f(find(ft>=t,1),d,d_); %= interp1(ft, f(:,d), t); % Interpolate the data set (ft, f) at times t
%         g0(d)=g(find(gt>=t,1),d);% = interp1(gt, g(:,d), t); % Interpolate the data set (ft, f) at times t
%        %f(d,d_)= interp1(ft, f(:,d,d_), t); % Interpolate the data set (ft, f) at times t
%        %g(d)=g(find(gt>=t,1),d);% = interp1(gt, g(:,d), t); % Interpolate the data set (ft, f) at times t
%     end
% end

f0=squeeze(f(find(ft>=t,1),:,:)); %= interp1(ft, f(:,d), t); % Interpolate the data set (ft, f) at times t
g0=squeeze(g(find(gt>=t,1),:));% = interp1(gt, g(:,d), t); % Interpolate the data set (ft, f) at times t
  
M=zeros(D,D);

M=f0;
for d=1:D
     M(d,d)=f0(d,d)-g0(d);  
end

dydt = -M*y; % Evalute ODE at times t
