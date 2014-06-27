%% mfileName
% (06/21/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)
%%
n = 500;
p = 200;
q = 2;

% noise level
SIGMA=5;

% noise correlation
rho=0.7;
%%  figure out how to sample from unif([-5,-3] || [3,5])
% K=500e3;
% xx= 2* rand(K,1) + 3;
% purge,figure,hist(xx,100)
% tmp_mask = randn(K,1)>0;
% xx(tmp_mask) = xx(tmp_mask)*-1;
% figure,hist(xx,100)
% 
% % one liner
% figure,hist(  (2*rand(K,1)+3) .* (  -1+2*round(rand(K,1))  ),       100)
%% create weight coefficients with random support
%=========================================================================%
% create support matrix
%=========================================================================%
% k=100;
k = round(0.25*p*q);
idx_supp = randsample(p*q,k);
mask=false(p,q);
mask(idx_supp)=true;
% figure,imagesc(mask),colormap(1-gray),imexpl,

%=========================================================================%
% create coefficients: unif([-5,-2} || [2,+5])
%=========================================================================%
Wtrue = zeros(p,q);
% Wtrue(mask) = rand(k,1)-0.5;
% Wtrue(mask) = randn(k,1);
bias=.5;mag=.5;Wtrue(mask) = (mag*rand(k,1)+bias) .*       (-1+2*round(rand(k,1)));
% figure,imexpt,subplot(121),imagesc(Wtrue),subplot(122),hist(Wtrue(mask),20),drawnow,return

%=========================================================================% 
% sample correlated noise
%=========================================================================%
COV = [1,rho;rho,1];
L=chol(COV)';

E = (L*randn(q,n))';

% figure,imexpb,tplot2(E(:,1),E(:,2))
% corr(E)
%=========================================================================%
% create data matrix
%=========================================================================%
X = randn(n,p);

%=========================================================================%
% sample data
%=========================================================================%
Y = X*Wtrue+SIGMA*E;
% figure,imexpb,tplot2(Y(:,1),Y(:,2)),return

%=========================================================================%
% test set
%=========================================================================%
ntest=100;
Xtest=randn(ntest,p);
Etest=(L*randn(q,ntest))';
Ytest=Xtest*Wtrue+SIGMA*Etest;

CORR_Y = corr(Y)
CORR_E = corr(E)
CORR_SIG = corr(X*Wtrue)
% corr(B)
% figure,imexpt,CAXIS=[-1,1];
% subplot(131),imcov(CORR_Y),caxis(CAXIS),title('Corr(Y)')
% subplot(132),imcov(CORR_E),caxis(CAXIS),title('Corr(NOISE)')
% subplot(133),imcov(CORR_SIG),caxis(CAXIS),title('Corr(signal)')
%% set options for admm
options.rho=1;
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 5e-4;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination conditio
%% regression
%=========================================================================%
% elastic-net
%=========================================================================%
lam1=1;
gam1=1;
[W.EN_STL,out.EN_STL]=tak_admm_EN_regr_STL(X,Y,lam1,gam1,options,Wtrue);

figure,imexp
subplot(241),tplot(log10(out.EN_STL.fval(2:end))), title('log10(function value)')
subplot(242),tplot(log10(out.EN_STL.wdist)), title('log10(||wtrue-west||)')
subplot(243),tplot(log10(out.EN_STL.rel_changevec)), title('log10(wnew-wold)')

subplot(245),tplot(Wtrue)
subplot(246),tplot(W.EN_STL)
subplot(247),tplot2(Wtrue,W.EN_STL), legend('wtrue','west'),title('')
subplot(248),tplot(abs(Wtrue-W.EN_STL)),title('|wtrue-west|')
%%
% purge
lwid=2;
figure,imexpt
subplot(211)
    stem(Ytest(:,1),'linewidth',lwid),xlim([1,length(Ytest(:,1))]), hold on
    stem(Xtest*W.EN_STL(:,1),'linewidth',lwid,'color','k')
    legend('true','EN-STL'), ylim([min(Ytest(:,1)),max(Ytest(:,1))]),grid on
subplot(212)
    stem(Ytest(:,2),'linewidth',lwid),xlim([1,length(Ytest(:,2))]), hold on
    stem(Xtest*W.EN_STL(:,2),'linewidth',lwid,'color','k')
    legend('true','EN-STL'), ylim([min(Ytest(:,1)),max(Ytest(:,2))]),grid on
drawnow

% mse_RR = norm(ytest - Xtest*w_RR)
% mse_EN = norm(ytest - Xtest*w_EN)
% mse_GN = norm(ytest - Xtest*w_GN)
% mse_FL = norm(ytest - Xtest*w_FL)
% corr_RR = corr(ytest, Xtest*w_RR)
CORR.EN_STL(1) = corr(Ytest(:,1), Xtest*W.EN_STL(:,1));
CORR.EN_STL(2) = corr(Ytest(:,2), Xtest*W.EN_STL(:,2))
