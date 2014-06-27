%% june19_2d_FL_GN
% (06/19/2014)
%=========================================================================%
% - 2-d fused lasso trial
%=========================================================================%
%%
clear all;
purge
randn('state',0)
n = 200;
nx=80;
ny=80;
p = nx*ny;

wtrue = zeros(nx,nx);
wtrue(21:40, 61:70) = 1;
wtrue(54:66, 21:33) = 3;
wtrue(11:25, 11:25) = 5; 
% imcovv(wtrue)
% return

sig = 2;
X = randn(n,p);
y = X*wtrue(:) + sig*randn(n,1);

ntest=200;
Xtest = randn(ntest,p);
ytest = Xtest*wtrue(:);
% tplot(y)
% tplot(X*wtrue(:))
% tplott(wtrue(:))
%% set options for admm
options.rho=1;
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 5e-4;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = true;    % <- track function values (may slow alg.)
%% estimate weight vector
%-------------------------------------------------------------------------%
% ridge regression
%-------------------------------------------------------------------------%
w_RR = tak_ridge_regression(X,y,.05);
%%
% purge
%-------------------------------------------------------------------------%
% enet
%-------------------------------------------------------------------------%
lam1 = 2;  % L1 penalty weight
gam1 = 10; % fused lasso penalty weight

[w_EN, output_EN]=tak_admm_EN_regr(X,y,lam1,gam1,options,wtrue(:));
%%
%-------------------------------------------------------------------------%
% fused lasso
%-------------------------------------------------------------------------%
lam2 = 1;  % L1 penalty weight
gam2 = 2; % fused lasso penalty weight
C=tak_diffmat_2d([nx,ny],0);
[w_FL,output_FL]=tak_admm_FL_regr_pcg(X,y,lam2,gam2,options,C,[],wtrue(:));
%% GraphNet using fista
lam3=1;
gam3=10;
options_FBS=options;
options_FBS.tol = 1e-3;
tau1=svds(X,1)^2;
% [tau2,~,flag]=eigs(C'*C,1);
tau2 = 8; % i know this upper-bound 
options_FBS.tau = 1/(tau1+gam3*tau2);
[w_GN,out_GN]=tak_GN_regr_FISTA(X,y,lam3,gam3,options_FBS,C,wtrue(:));

figure,imexp
subplot(241),tplot(log10(out_GN.fval(2:end))), title('log10(function value)')
subplot(242),tplot(log10(out_GN.wdist)), title('log10(||wtrue-west||)')
subplot(243),tplot(log10(out_GN.rel_changevec)), title('log10(wnew-wold)')

subplot(245),tplot(wtrue(:))
subplot(246),tplot(w_GN(:))
subplot(247),tplot2(wtrue(:),w_GN(:)), legend('wtrue','west'),title('')
subplot(248),tplot(abs(wtrue(:)-w_GN(:))),title('|wtrue-west|')
% return
%%% compare methods
% purge
w_RR_array = reshape(w_RR, [nx,ny]);
w_EN_array = reshape(w_EN, [nx,ny]);
w_FL_array = reshape(w_FL, [nx,ny]);
w_GN_array = reshape(w_GN, [nx,ny]);

% figure,imexp
% subplot(221),imcov(wtrue)
% subplot(222),imcov(w_RR_array)
% subplot(223),imcov(w_EN_array)
% subplot(224),imcov(w_FL_array)

CAXIS=[-5,5];
figure,imexpt
% subplot(221),imcov(wtrue),caxis(CAXIS); colorbar('location','northoutside')
subplot(141),imcov(wtrue),CAXIS=caxis; colorbar('location','northoutside')
subplot(142),imcov(w_EN_array),caxis(CAXIS); colorbar('location','northoutside')
subplot(143),imcov(w_GN_array),caxis(CAXIS); colorbar('location','northoutside')
subplot(144),imcov(w_FL_array),caxis(CAXIS); colorbar('location','northoutside')
drawnow
        
lwid=2;
figure,imexpt
% subplot(221),imcov(wtrue),caxis(CAXIS); colorbar('location','northoutside')
stem(ytest,'linewidth',lwid),xlim([1,length(ytest)]), hold on
stem(Xtest*w_EN,'linewidth',lwid,'color','k')
stem(Xtest*w_GN,'linewidth',lwid,'color','g')
stem(Xtest*w_FL,'linewidth',lwid,'color','r')
legend('true','EN','GN','FL'), ylim([min(ytest),max(ytest)]),grid on
drawnow

% mse_RR = norm(ytest - Xtest*w_RR)
% mse_EN = norm(ytest - Xtest*w_EN)
% mse_GN = norm(ytest - Xtest*w_GN)
% mse_FL = norm(ytest - Xtest*w_FL)
% corr_RR = corr(ytest, Xtest*w_RR)
corr_EN = corr(ytest, Xtest*w_EN)
corr_GN = corr(ytest, Xtest*w_GN)
corr_FL = corr(ytest, Xtest*w_FL)