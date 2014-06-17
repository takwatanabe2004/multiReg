%% june14_1d_flasso
% (06/14/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
% purge
randn('state',0)
n = 888;
p = 3000;

wtrue = zeros(p,1);
wtrue(21:40) = 1;
wtrue(81:100) = 3;
wtrue(101:120) = 5; 
wtrue(121:130) = 2;
wtrue(221:230) = 2;
% data = rand(p,1); data(data>0.8) = 1; data(data<=0.8)=0;
% data = binornd(1,0.2,[p,1]);
% data(data==1)

sig = 3;
X = randn(n,p);
y = X*wtrue + sig*randn(n,1);
% tplot(y)
% tplot(X*wtrue)
% tplott(wtrue)

%% ridge regression
west = tak_ridge_regression(X,y,.05);
%% enet
options.lambda=99;  % L1 penalty weight
options.gamma =0; % fused lasso penalty weight
options.rho=1;
% termination criterion
options.termin.maxiter = 500;   % <- maximum number of iterations
options.termin.tol = 5e-4;      % <- relative change in the primal variable
options.termin.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.termin.silence = false; % <- display termination condition

west_L1=tak_admm_enet_regr(X,y,options,wtrue);
% cppd
clear options
options.lambda= 55;  % L1 penalty weight
options.gamma = 0; % fused lasso penalty weight

% termination criterion
options.termin.maxiter = 500;   % <- maximum number of iterations
options.termin.tol = 5e-4;      % <- relative change in the primal variable
options.termin.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.termin.silence = false; % <- display termination condition
options.fval = true;
C = tak_diffmat_1d(p,0);

options.sigma=1;
F = [options.lambda*speye(p);options.gamma*C];
% tic
% options.L=svds(F,1)
% options.L=sqrt(eigs(F'*F,1))
options.L = normest(F'*F)
% toc
options.tau=1/(options.L^2 * options.sigma)
options.tau = options.tau - options.tau/100; % <- safeguard (sig*tau*L^2 < 1...strict equality)
% return
% output=tak_cppd_flas_regr(X,y,options,C,wtrue)
west_cppd=output.w;

purge
figure,imexp
subplot(311),tstairs2(wtrue,west),ylim([-1,6.5])
subplot(312),tstairs2(wtrue,west_L1.w),ylim([-1,6.5])
subplot(313),tstairs2(wtrue,west_cppd),ylim([-1,6.5])
% isequal(west, tak_ridge_regression(X,y))
drawnow

