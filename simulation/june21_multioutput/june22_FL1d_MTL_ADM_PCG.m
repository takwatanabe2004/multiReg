%% june22_FL1d_MTL_ADM_PCG
% (06/22/2014)
%=========================================================================%
% - ADMM-PCG implementation of multitask Fused Lasso
% - ground truth has 1-d patch structure
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)
%%
n = 200;
p = 300;
q=50;

% noise level
SIGMA=111;

% penalty values
lam_stl=1533;
gam_stl=25;
lam_mtl=6811;
gam_mtl=gam_stl;

% for FL
lamFL_stl=lam_stl;
gamFL_stl=55;
lamFL_mtl=lam_mtl;
gamFL_mtl=gamFL_stl;
%%
%=========================================================================%
% create weight vector with intra-piecewise-constant-1d structure, that has shared
% sparsity pattern across outputs (row sparsity)
%=========================================================================%
Wtrue = zeros(p,q);
k=25;
for col=1:q
    % random amplitude
    amp = (3+20*rand(5,1)).*sign(randn(5,1));    
    idx1 = (1:k)+2*k;
    idx2 = (1:k)+3*k;
    idx3 = (1:k)+4*k;
    idx4 = (1:k)+p-k;
    idx5 = (1:k)+p-2*k;
%     Wtrue(1+k:2*k,col)   = amp(2);
    Wtrue(idx1,col)   = amp(1);
    Wtrue(idx2,col)   = amp(2);
    Wtrue(idx3,col)   = amp(3);
    Wtrue(idx4,col)   = amp(4);
    Wtrue(idx5,col)   = amp(5);
%     Wtrue(501:500+k,col) = amp(3);
%     Wtrue(601:600+k,col) = amp(4);
end
% figure,imexpl,plot(Wtrue(:,1:3),'linewidth',3),grid on,
% figure,imexpl,imagesc(Wtrue)
% figure,imexpl,imagesc(Wtrue~=0)
% drawnow,return

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
% figure,imexpb,plot(Y(:,1:5),'linewidth',2),grid on,drawnow,return

%=========================================================================%
% test set
%=========================================================================%
ntest=100;
Xtest=randn(ntest,p);
Etest=randn(ntest,q);
Ytest=Xtest*Wtrue+SIGMA*Etest;
%% set options
options.rho=.5;
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 5e-6;      % <- relative change in the primal variable
options.progress = 100;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = true;    % <- track function values (may slow alg.)

C=tak_diffmat_1d(p,0);

%-------------------------------------------------------------------------%
% spectral radiuses (needed for fista stepsize)
%-------------------------------------------------------------------------%
tau_X=svds(X,1)^2;
tau_C = tnormest(C'*C);
% return
%% regression
%=========================================================================%
% EN-STL: ADM
%=========================================================================%
options.MTL=false;
[W.EN_STL,out.EN_STL]=...
    tak_EN_regr_MTL_ADMM(X,Y,lam_stl,gam_stl,options,Wtrue);

figure,imexp
subplot(231),tplot(log10(out.EN_STL.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.EN_STL.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.EN_STL.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.EN_STL),mytitle('EN_STL'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.EN_STL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
% return
%=========================================================================%
% EN-MTL: ADM
%=========================================================================%
options.MTL=true;
[W.EN_MTL,out.EN_MTL]=...
    tak_EN_regr_MTL_ADMM(X,Y,lam_mtl,gam_mtl,options,Wtrue);

figure,imexp
subplot(231),tplot(log10(out.EN_MTL.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.EN_MTL.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.EN_MTL.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.EN_MTL),mytitle('EN_MTL'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.EN_MTL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
% % return

%=========================================================================%
% FL-STL: ADM
%=========================================================================%
options.MTL=false;
[W.FL_STL,out.FL_STL]=...
    tak_FL_regr_MTL_ADMM_PCG(X,Y,lamFL_stl,gamFL_stl,options,C,[],Wtrue);

figure,imexp
subplot(231),tplot(log10(out.FL_STL.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.FL_STL.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.FL_STL.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.FL_STL),mytitle('FL_STL'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.FL_STL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
% return
%=========================================================================%
% FL-MTL: ADM
%=========================================================================%
options.MTL=true;
[W.FL_MTL,out.FL_MTL]=...
    tak_FL_regr_MTL_ADMM_PCG(X,Y,lamFL_mtl,gamFL_mtl,options,C,[],Wtrue);

figure,imexp
subplot(231),tplot(log10(out.FL_MTL.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.FL_MTL.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.FL_MTL.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.FL_MTL),mytitle('FL_MTL'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.FL_MTL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
% return
%% regression: EN-fista, GN-fista
% lam_stl=33;
% gam_stl=5;
% lam_mtl=222;
% gam_mtl=gam_stl;
% 
% %=========================================================================%
% % EN-STL
% %=========================================================================%
% % get stepsize for ista
% options.tau = 1/(tau_X+gam_stl);
% 
% options.MTL=false;
% [W.EN_STL,out.EN_STL]=...
%     tak_EN_regr_MTL_FISTA(X,Y,lam_stl,gam_stl,options,Wtrue);
% 
% figure,imexp
% subplot(231),tplot(log10(out.EN_STL.fval(2:end))), title('log10(function value)')
% subplot(232),tplot(log10(out.EN_STL.wdist)), title('log10(||wtrue-west||)')
% subplot(233),tplot(log10(out.EN_STL.rel_changevec)), title('log10(wnew-wold)')
% subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(235),imagesc(W.EN_STL),mytitle('EN_STL'),caxis(CAXIS),colorbar
% subplot(236),imagesc(abs(Wtrue-W.EN_STL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
% drawnow
% % return
% 
% %=========================================================================%
% % EN-MTL:
% %=========================================================================%
% options.MTL=true;
% options.tau = 1/(tau_X+gam_mtl);
% [W.EN_MTL,out.EN_MTL]=...
%     tak_EN_regr_MTL_FISTA(X,Y,lam_mtl,gam_mtl,options,Wtrue);
% 
% figure,imexp
% subplot(231),tplot(log10(out.EN_MTL.fval(2:end))), title('log10(function value)')
% subplot(232),tplot(log10(out.EN_MTL.wdist)), title('log10(||wtrue-west||)')
% subplot(233),tplot(log10(out.EN_MTL.rel_changevec)), title('log10(wnew-wold)')
% subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(235),imagesc(W.EN_MTL),mytitle('EN_MTL'),caxis(CAXIS),colorbar
% subplot(236),imagesc(abs(Wtrue-W.EN_MTL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
% drawnow
% % return
% 
%=========================================================================%
% GN-STL
%=========================================================================%
options.MTL=false;
options.tau = 1/(tau_X+gam_stl*tau_C);
[W.GN_STL,out.GN_STL]=...
    tak_GN_regr_MTL_FISTA(X,Y,lam_stl,gam_stl,options,C,Wtrue);
figure,imexp
subplot(231),tplot(log10(out.GN_STL.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.GN_STL.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.GN_STL.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.GN_STL),mytitle('GN_STL'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.GN_STL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
% return

%=========================================================================%
% GN-MTL
%=========================================================================%
options.MTL=true;
options.tau = 1/(tau_X+gam_mtl*tau_C);
[W.GN_MTL,out.GN_MTL]=...
    tak_GN_regr_MTL_FISTA(X,Y,gam_mtl,gam_mtl,options,C,Wtrue);
figure,imexp
subplot(231),tplot(log10(out.GN_MTL.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.GN_MTL.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.GN_MTL.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.GN_MTL),mytitle('GN_MTL'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.GN_MTL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
% % return
%% evaluate test performance
Ypred.EN_STL = Xtest*W.EN_STL;
Ypred.EN_MTL = Xtest*W.EN_MTL;
Ypred.GN_STL = Xtest*W.GN_STL;
Ypred.GN_MTL = Xtest*W.GN_MTL;
Ypred.FL_STL = Xtest*W.FL_STL;
Ypred.FL_MTL = Xtest*W.FL_MTL;

MSE.EN_STL = norm(Ytest(:) - Ypred.EN_STL(:));
MSE.EN_MTL = norm(Ytest(:) - Ypred.EN_MTL(:));
MSE.GN_STL = norm(Ytest(:) - Ypred.GN_STL(:));
MSE.GN_MTL = norm(Ytest(:) - Ypred.GN_MTL(:));
MSE.FL_STL = norm(Ytest(:) - Ypred.FL_STL(:));
MSE.FL_MTL = norm(Ytest(:) - Ypred.FL_MTL(:));

CORR.EN_STL=corr(Ytest(:), Ypred.EN_STL(:));
CORR.EN_MTL=corr(Ytest(:), Ypred.EN_MTL(:));
CORR.GN_STL=corr(Ytest(:), Ypred.GN_STL(:));
CORR.GN_MTL=corr(Ytest(:), Ypred.GN_MTL(:));
CORR.FL_STL=corr(Ytest(:), Ypred.FL_STL(:));
CORR.FL_MTL=corr(Ytest(:), Ypred.FL_MTL(:));

MSE
CORR
% 
% drawnow