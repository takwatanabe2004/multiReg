%% june16_2d_flasso
% (06/16/2014)
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
options.funcval = false; 
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
lam1 = 25;  % L1 penalty weight
gam1 = 10; % fused lasso penalty weight

[w_EN, output_EN]=tak_admm_EN_regr(X,y,lam1,gam1,options,wtrue(:));
%%
%-------------------------------------------------------------------------%
% fused lasso
%-------------------------------------------------------------------------%
lam2 = .5;  % L1 penalty weight
gam2 = 1; % fused lasso penalty weight
C=tak_diffmat_2d([nx,ny],0);
[w_FL,output_FL]=tak_admm_FL_regr_pcg(X,y,lam2,gam2,options,C,[],wtrue(:));

% figure,imexp
% subplot(241),tplot(log10(output.fval(2:end))), title('log10(function value)')
% subplot(242),tplot(log10(output.wdist)), title('log10(||wtrue-west||)')
% subplot(243),tplot(log10(output.rel_changevec)), title('log10(wnew-wold)')
% 
% subplot(245),tplot(wtrue)
% subplot(246),tplot(w_FL)
% subplot(247),tplot2(wtrue,w_FL), legend('wtrue','west'),title('')
% subplot(248),tplot(abs(wtrue-w_FL)),title('|wtrue-west|')

%%% compare methods
% purge
w_RR_array = reshape(w_RR, [nx,ny]);
w_EN_array = reshape(w_EN, [nx,ny]);
w_FL_array = reshape(w_FL, [nx,ny]);

% figure,imexp
% subplot(221),imcov(wtrue)
% subplot(222),imcov(w_RR_array)
% subplot(223),imcov(w_EN_array)
% subplot(224),imcov(w_FL_array)

CAXIS=[-5,5];
figure,imexpb
% subplot(221),imcov(wtrue),caxis(CAXIS); colorbar('location','northoutside')
subplot(141),imcov(wtrue),CAXIS=caxis; colorbar('location','northoutside')
subplot(142),imcov(w_RR_array),caxis(CAXIS); colorbar('location','northoutside')
subplot(143),imcov(w_EN_array),caxis(CAXIS); colorbar('location','northoutside')
subplot(144),imcov(w_FL_array),caxis(CAXIS); colorbar('location','northoutside')
drawnow
        
lwid=2;
figure,imexpb
% subplot(221),imcov(wtrue),caxis(CAXIS); colorbar('location','northoutside')
stem(ytest,'linewidth',lwid),xlim([1,length(ytest)]), hold on
stem(Xtest*w_RR,'linewidth',lwid,'color','k')
stem(Xtest*w_EN,'linewidth',lwid,'color','g')
stem(Xtest*w_FL,'linewidth',lwid,'color','r')
legend('true','RR','EN','FL'), ylim([min(ytest),max(ytest)]),grid on
drawnow

mse_RR = norm(ytest - Xtest*w_RR)
mse_EN = norm(ytest - Xtest*w_EN)
mse_FL = norm(ytest - Xtest*w_FL)
corr_RR = corr(ytest, Xtest*w_RR)
corr_EN = corr(ytest, Xtest*w_EN)
corr_FL = corr(ytest, Xtest*w_FL)