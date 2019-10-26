function dydt = CVM_CTBN_mu_fastD(t, y,D, ft, f,gt,g)

%%% AS ODE IS TIME DEPENDENT ON EFF. RATES WE HAVE TO INTERPOLATE



f0=squeeze(f(find(ft>=t,1),:,:)); %= interp1(ft, f(:,d), t); % Interpolate the data set (ft, f) at times t
g0=squeeze(g(find(gt>=t,1),:));% = interp1(gt, g(:,d), t); % Interpolate the data set (ft, f) at times t
       % f(d,d_)=interp1(ft, f(:,d,d_), t); % Interpolate the data set (ft, f) at times t
       % g(d)=g(find(gt>=t,1),d);% = interp1(gt, g(:,d), t); 

M=zeros(D,D);

% for d=1:D
%     for d_=1:D
%         if d~=d_
%             M(d,d_)=f(d_,d)*g(d)/g(d_);
%         end
%     end
% end

M=f0'.*g0'./g0;
%M=f'*g/g';

for d=1:D
    M(d,d)=0;  
    M(d,d)=-sum(M(:,d));
end

%M=[-f(1)*g(2)/g(1) +f(1)*g(2)/g(1) ; f(1)*g(2)/g(1) -f(2)*g(1)/g(2)];
%M=[-f(1)*g(1)/g(2) +f(2)*g(2)/g(1) ; f(1)*g(1)/g(2) -f(2)*g(2)/g(1)];
dydt = M*y; % Evalute ODE at times t
