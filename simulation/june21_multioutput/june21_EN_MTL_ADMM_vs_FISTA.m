%% june21_MTL_ADMM_vs_FISTA
% (06/21/2014)
%=========================================================================%
% - compare EN implementation: ADMM vs FISTA
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)
%%
n = 500;
p = 600;
q = 80;

% noise level
SIGMA=2;
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
options.rho=45;
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 1e-4;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = true;    % <- track function values (may slow alg.)
%% regression
%=========================================================================%
% elastic-net-ista
%=========================================================================%
lam2=333;
gam2=0;

% get stepsize for ista
tau=svds(X,1)^2;
options.tau = 1/(tau+gam2);

% run ista
[W.EN_IST,out.EN_IST]=tak_EN_regr_MTL_ISTA(X,Y,lam2,gam2,options,Wtrue);

% figure,imexp
% subplot(231),tplot(log10(out.EN_IST.fval(2:end))), title('log10(function value)')
% subplot(232),tplot(log10(out.EN_IST.wdist)), title('log10(||wtrue-west||)')
% subplot(233),tplot(log10(out.EN_IST.rel_changevec)), title('log10(wnew-wold)')
% subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(235),imagesc(W.EN_IST),title('West-EN-ISTA'),caxis(CAXIS),colorbar
% subplot(236),imagesc(abs(Wtrue-W.EN_IST)),title('|wtrue-west|'),caxis(CAXIS),colorbar

%=========================================================================%
% elastic-net-fista
%=========================================================================%
[W.EN_FIST,out.EN_FIST]=tak_EN_regr_MTL_FISTA(X,Y,lam2,gam2,options,Wtrue);

% figure,imexp
% subplot(231),tplot(log10(out.EN_FIST.fval(2:end))), title('log10(function value)')
% subplot(232),tplot(log10(out.EN_FIST.wdist)), title('log10(||wtrue-west||)')
% subplot(233),tplot(log10(out.EN_FIST.rel_changevec)), title('log10(wnew-wold)')
% subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(235),imagesc(W.EN_FIST),title('West-EN-FISTA'),caxis(CAXIS),colorbar
% subplot(236),imagesc(abs(Wtrue-W.EN_FIST)),title('|wtrue-west|'),caxis(CAXIS),colorbar

%=========================================================================%
% elastic-net-admm
%=========================================================================%
lam2=333;
gam2=0;
[W.EN_ADM,out.EN_ADM]=tak_EN_regr_MTL_ADMM(X,Y,lam2,gam2,options,Wtrue);

% figure,imexp
% subplot(231),tplot(log10(out.EN_ADM.fval(2:end))), title('log10(function value)')
% subplot(232),tplot(log10(out.EN_ADM.wdist)), title('log10(||wtrue-west||)')
% subplot(233),tplot(log10(out.EN_ADM.rel_changevec)), title('log10(wnew-wold)')
% subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(235),imagesc(W.EN_ADM),title('West-EN-ADMM'),caxis(CAXIS),colorbar
% subplot(236),imagesc(abs(Wtrue-W.EN_ADM)),title('|wtrue-west|'),caxis(CAXIS),colorbar
%% compare
% purge
figure,imexp
subplot(241)
    tplot3(log10(out.EN_ADM.fval),  log10(out.EN_IST.fval),...
           log10(out.EN_FIST.fval))
    HH=legend('ADMM','ISTA','FISTA'); set(HH,'Interpreter','none')
    title('log10(function value)');
subplot(242),
    tplot3(log10(out.EN_ADM.wdist),  log10(out.EN_IST.wdist),...
           log10(out.EN_FIST.wdist))
    HH=legend('ADMM','ISTA','FISTA'); set(HH,'Interpreter','none')
    title('log10(||wtrue-west||)')
subplot(243),
    tplot3(log10(out.EN_ADM.rel_changevec),  log10(out.EN_IST.rel_changevec),...
           log10(out.EN_FIST.rel_changevec))
    HH=legend('ADMM','ISTA','FISTA'); set(HH,'Interpreter','none')
    title('log10(wnew-wold)')
% subplot(244),tplot3(w_EN_ADM,w_EN_FIST_ver2,w_EN_FIST), 
%     HH=legend('w_est_ADM','w_GN_FIST_ver2','w_est_FIST');set(HH,'Interpreter','none')
subplot(245),imagesc(Wtrue)
subplot(245),imagesc(W.EN_ADM),title('ADM')
subplot(246),imagesc(W.EN_IST),title('ISTA')
subplot(247),imagesc(W.EN_FIST),title('FISTA')
subplot(248),
    imagesc(abs(Wtrue-W.EN_FIST)),title('|wtrue-w_FISTA|','Interpreter','none')

figure,imexp
subplot(131)
    tplot2(log10(out.EN_IST.fval),...
           log10(out.EN_FIST.fval))
    HH=legend('ISTA','FISTA'); set(HH,'Interpreter','none')
    title('log10(function value)');
subplot(132),
    tplot2(log10(out.EN_IST.wdist),...
           log10(out.EN_FIST.wdist))
    HH=legend('ISTA','FISTA'); set(HH,'Interpreter','none')
    title('log10(||wtrue-west||)')
subplot(133),
    tplot2(log10(out.EN_IST.rel_changevec),...
           log10(out.EN_FIST.rel_changevec))
    HH=legend('ISTA','FISTA'); set(HH,'Interpreter','none')
    title('log10(wnew-wold)')
%%
% Ypred.RR = Xtest*W.RR;
% Ypred.EN_STL = Xtest*W.EN_STL;
% Ypred.EN_ADM = Xtest*W.EN_ADM;
% 
% MSE.RR = norm(Ytest(:) - Ypred.RR(:));
% MSE.EN_STL = norm(Ytest(:) - Ypred.EN_STL(:));
% MSE.EN_ADM = norm(Ytest(:) - Ypred.EN_ADM(:));
% CORR.RR=corr(Ytest(:), Ypred.RR(:));
% CORR.EN_STL=corr(Ytest(:), Ypred.EN_STL(:));
% CORR.EN_ADM=corr(Ytest(:), Ypred.EN_ADM(:));
% 
% MSE
% CORR
% 
% drawnow