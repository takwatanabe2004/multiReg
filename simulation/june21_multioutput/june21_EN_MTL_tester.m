%% june21_MTL_tester
% (06/21/2014)
%=========================================================================%
% - here common support is shared among the weight coefficients
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)
%%
n = 100;
p = 60;
q = 80;

% noise level
SIGMA=5;
%% create weight coefficients with random support
%=========================================================================%
% create support matrix with shared support (row sparsity)
%=========================================================================%
% k=100;
k = round(0.1*p);
idx_supp = randsample(p,k);
mask=false(p,q);
mask(idx_supp,:)=true;
% figure,imagesc(mask),colormap(1-gray),imexpl,drawnow,return

%=========================================================================%
% create coefficients: unif([-5,-2} || [2,+5])
%=========================================================================%
Wtrue = zeros(p,q);
% Wtrue(mask) = rand(k,q)-0.5;
% Wtrue(mask) = randn(k,q);
bias=.5;mag=.5;Wtrue(mask) = (mag*rand(k,q)+bias) .*       (-1+2*round(rand(k,q)));
% figure,imexpt,subplot(121),imagesc(Wtrue),subplot(122),hist(Wtrue(mask),20),drawnow,return

%=========================================================================% 
% sample white noise
%=========================================================================%
E = randn(n,q);

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
Etest=randn(ntest,q);
Ytest=Xtest*Wtrue+SIGMA*Etest;
%% set options for admm
options.rho=1;
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 1e-3;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination conditio
%% regression
%=========================================================================%
% elastic-net-stl
%=========================================================================%
lam1=22;
gam1=0;
[W.EN_STL,out.EN_STL]=tak_admm_EN_regr_STL(X,Y,lam1,gam1,options,Wtrue);

figure,imexp
subplot(231),tplot(log10(out.EN_STL.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.EN_STL.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.EN_STL.rel_changevec)), title('log10(wnew-wold)')

subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.EN_STL),title('West-EN-STL'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.EN_STL)),title('|wtrue-west|'),caxis(CAXIS),colorbar

%=========================================================================%
% elastic-net-mtl
%=========================================================================%
lam2=333;
gam2=0;
[W.EN_MTL,out.EN_MTL]=tak_admm_EN_regr_MTL(X,Y,lam2,gam2,options,Wtrue);

figure,imexp
subplot(231),tplot(log10(out.EN_MTL.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.EN_MTL.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.EN_MTL.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.EN_MTL),title('West-EN-MTL'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.EN_MTL)),title('|wtrue-west|'),caxis(CAXIS),colorbar

%=========================================================================%
% ridge 
%=========================================================================%
W.RR = tak_ridge_regression(X,Y,lam1);
% return
figure,imexpt
subplot(141),imagesc(Wtrue),caxis(CAXIS),colorbar,title('Wtrue')
subplot(142),imagesc(W.RR),caxis(CAXIS),colorbar,title('RR')
subplot(143),imagesc(W.EN_STL),caxis(CAXIS),colorbar,title('EN-STL')
subplot(144),imagesc(W.EN_MTL),caxis(CAXIS),colorbar,title('EN-MTL')
%%
Ypred.RR = Xtest*W.RR;
Ypred.EN_STL = Xtest*W.EN_STL;
Ypred.EN_MTL = Xtest*W.EN_MTL;

MSE.RR = norm(Ytest(:) - Ypred.RR(:));
MSE.EN_STL = norm(Ytest(:) - Ypred.EN_STL(:));
MSE.EN_MTL = norm(Ytest(:) - Ypred.EN_MTL(:));
CORR.RR=corr(Ytest(:), Ypred.RR(:));
CORR.EN_STL=corr(Ytest(:), Ypred.EN_STL(:));
CORR.EN_MTL=corr(Ytest(:), Ypred.EN_MTL(:));

MSE
CORR

drawnow