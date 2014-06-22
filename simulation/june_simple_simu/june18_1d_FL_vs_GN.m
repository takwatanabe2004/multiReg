%% june18_1d_FL_vs_GN
% (06/18/2014)
%=========================================================================%
% - My protocode for GN-penalized regression
%=========================================================================%
%%
clear all;
% purge
randn('state',0)
n = 100;
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

sig = .2;
X = randn(n,p);
y = X*wtrue + sig*randn(n,1);
% tplot(y)
% tplot(X*wtrue)
% tplott(wtrue)
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
lam2 = 1;  % L1 penalty weight
gam2 = 2; % fused lasso penalty weight
[w_FL,output_FL]=tak_admm_FL_regr_pcg(X,y,lam2,gam2,options,tak_diffmat_1d(p,0),[],wtrue);

%-------------------------------------------------------------------------%
% graphnet
%-------------------------------------------------------------------------%
lam3 = 0;  % L1 penalty weight
gam3 = 1; % fused lasso penalty weight
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
subplot(411),tstairs2(wtrue,w_RR),ylim([-1,6.5])
subplot(412),tstairs2(wtrue,w_EN),ylim([-1,6.5])
subplot(413),tstairs2(wtrue,w_FL),ylim([-1,6.5])
subplot(414),tstairs2(wtrue,w_GN),ylim([-1,6.5])
drawnow
