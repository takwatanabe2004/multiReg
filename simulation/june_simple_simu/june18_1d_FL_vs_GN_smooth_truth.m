%% june18_1d_FL_vs_GN_smooth_truth
% (06/18/2014)
%=========================================================================%
% - same as june18_1d_FL_vs_GN.m, but with smooth weight vector
%=========================================================================%
%%
clear all;
% purge
randn('state',0)
n = 100;
p = 2000;

k=200;
wtrue = zeros(p,1);
wtrue(1:k) = 4*sin(.02*pi*(1:k));
wtrue(601:600+k) = 4*sin(.02*pi*(1:k));
% data = rand(p,1); data(data>0.8) = 1; data(data<=0.8)=0;
% data = binornd(1,0.2,[p,1]);
% data(data==1)

sig = .2;
X = randn(n,p);
y = X*wtrue + sig*randn(n,1);

ntest=200;
Xtest=randn(ntest,p);
ytest=Xtest*wtrue;

% tplot(y)
% tplot(X*wtrue)
% tplott(wtrue)
% return
%% set options for admm
options.rho=1;
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 5e-4;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = false;    % <- track function values (may slow alg.)
%% estimate weight vector
%-------------------------------------------------------------------------%
% ridge regression
%-------------------------------------------------------------------------%
w_RR = tak_ridge_regression(X,y,.05);

%-------------------------------------------------------------------------%
% enet
%-------------------------------------------------------------------------%
lam1 = 3;  % L1 penalty weight
gam1 = 0; % fused lasso penalty weight

[w_EN, output_EN]=tak_admm_EN_regr(X,y,lam1,gam1,options,wtrue);

%-------------------------------------------------------------------------%
% fused lasso
%-------------------------------------------------------------------------%
lam2 = .25;  % L1 penalty weight
gam2 = 8; % fused lasso penalty weight
[w_FL,output_FL]=tak_admm_FL_regr_pcg(X,y,lam2,gam2,options,tak_diffmat_1d(p,0),[],wtrue);

%-------------------------------------------------------------------------%
% graphnet
%-------------------------------------------------------------------------%
lam3 = 1;  % L1 penalty weight
gam3 = 23; % fused lasso penalty weight
[w_GN,output_GN]=tak_admm_GN_regr_pcg(X,y,lam3,gam3,options,tak_diffmat_1d(p,0),[],wtrue);

% figure,imexp
% subplot(241),tplot(log10(output.fval(2:end))), title('log10(function value)')
% subplot(242),tplot(log10(output.wdist)), title('log10(||wtrue-west||)')
% subplot(243),tplot(log10(output.rel_changevec)), title('log10(wnew-wold)')
% 
% subplot(245),tplot(wtrue)
% subplot(246),tplot(w_FL)
% subplot(247),tplot2(wtrue,w_FL), legend('wtrue','west'),title('')
% subplot(248),tplot(abs(wtrue-w_FL)),title('|wtrue-west|')
%% compare methods
% purge
figure,imexp
subplot(411),stairs(wtrue),YLIM=ylim; YLIM = [YLIM(1)-1,YLIM(2)+1];
subplot(411),tstairs2(wtrue,w_RR),ylim(YLIM)
subplot(412),tstairs2(wtrue,w_EN),ylim(YLIM)
subplot(413),tstairs2(wtrue,w_FL),ylim(YLIM)
subplot(414),tstairs2(wtrue,w_GN),ylim(YLIM)
drawnow

mse_RR = norm(ytest - Xtest*w_RR)
mse_EN = norm(ytest - Xtest*w_EN)
mse_FL = norm(ytest - Xtest*w_FL)
mse_GN = norm(ytest - Xtest*w_GN)

corr_RR = corr(ytest , Xtest*w_RR)
corr_EN = corr(ytest , Xtest*w_EN)
corr_FL = corr(ytest , Xtest*w_FL)
corr_GN = corr(ytest , Xtest*w_GN)