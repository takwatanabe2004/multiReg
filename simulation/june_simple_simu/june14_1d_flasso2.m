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
lam=99;  % L1 penalty weight
gam =0; % fused lasso penalty weight
options.rho=1;
% termination criterion
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 5e-4;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition

west_EN=tak_admm_EN_regr(X,y,lam,gam,[],wtrue);

purge
figure,imexp
subplot(311),tstairs2(wtrue,west),ylim([-1,6.5])
subplot(312),tstairs2(wtrue,west_EN),ylim([-1,6.5])
% subplot(313),tstairs2(wtrue,west_cppd),ylim([-1,6.5])
% isequal(west, tak_ridge_regression(X,y))
drawnow

