%% june16_1d_flasso
% (06/16/2014)
%=========================================================================%
% - 1-d fused lasso trial
%=========================================================================%
%%
clear all;
% purge
randn('state',0)
n = 50;
p = 1000;

wtrue = zeros(p,1);
wtrue(21:40) = 1;
wtrue(81:100) = 3;
wtrue(101:120) = 5; 
wtrue(121:130) = 2;
wtrue(221:230) = 2;
% data = rand(p,1); data(data>0.8) = 1; data(data<=0.8)=0;
% data = binornd(1,0.2,[p,1]);
% data(data==1)

sig = 1;
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
%% fused lasso
lam = 1;  % L1 penalty weight
gam = 2; % fused lasso penalty weight

options.rho=.1;
% termination criterion
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 1e-4;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition

C = tak_diffmat_1d(p,0);
% output=tak_admm_FL_regr_pcg(X,y,lam,gam,options,C,PCG,wtrue)
output=tak_admm_FL_regr_pcg(X,y,lam,gam,options,C,[],wtrue)
west_fl=output.w;



% figure,imexp
% subplot(241),tplot(log10(output.fval(2:end))), title('log10(function value)')
% subplot(242),tplot(log10(output.wdist)), title('log10(||wtrue-west||)')
% subplot(243),tplot(log10(output.rel_changevec)), title('log10(wnew-wold)')
% 
% subplot(245),tplot(wtrue)
% subplot(246),tplot(west_fl)
% subplot(247),tplot2(wtrue,west_fl), legend('wtrue','west'),title('')
% subplot(248),tplot(abs(wtrue-west_fl)),title('|wtrue-west|')
%% compare methods
% purge
figure,imexp
subplot(311),tstairs2(wtrue,west),ylim([-1,6.5])
subplot(312),tstairs2(wtrue,west_L1.w),ylim([-1,6.5])
subplot(313),tstairs2(wtrue,west_fl),ylim([-1,6.5])
drawnow
