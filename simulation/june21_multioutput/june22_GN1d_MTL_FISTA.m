%% june22_GN1d_MTL_FISTA
% (06/22/2014)
%=========================================================================%
% - FISTA implementation of multitask GraphNet
% - ground truth has 1-d smooth structure
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)
%%
n = 200;
p = 1000;
q=200;

% noise level
SIGMA=20;
%%
%=========================================================================%
% create weight vector with intra-smooth-1d structure, that has shared
% sparsity pattern across outputs (row sparsity)
%=========================================================================%
Wtrue = zeros(p,q);
k=50;
for col=1:q
    % random amplitude
    amp = (2+5*rand(4,1)).*sign(randn(4,1));    
    Wtrue(1:k,col)       = amp(1)*sin(.02*pi*(1:k));    
    Wtrue(1+k:2*k,col)   = amp(2)*sin(.02*pi*(1:k));    
    Wtrue(501:500+k,col) = amp(3)*sin(.02*pi*(1:k));    
    Wtrue(601:600+k,col) = amp(4)*sin(.02*pi*(1:k));
end
% figure,imexpl,plot(Wtrue(:,1:5),'linewidth',3),grid on,
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
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 1e-3;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
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
lam_stl=33;
gam_stl=5;
lam_mtl=222;
gam_mtl=gam_stl;

%=========================================================================%
% EN-STL
%=========================================================================%
% get stepsize for ista
options.tau = 1/(tau_X+gam_stl);

options.MTL=false;
[W.EN_STL,out.EN_STL]=...
    tak_EN_regr_MTL_FISTA(X,Y,lam_stl,gam_stl,options,Wtrue);

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
% EN-MTL:
%=========================================================================%
options.MTL=true;
options.tau = 1/(tau_X+gam_mtl);
[W.EN_MTL,out.EN_MTL]=...
    tak_EN_regr_MTL_FISTA(X,Y,lam_mtl,gam_mtl,options,Wtrue);

figure,imexp
subplot(231),tplot(log10(out.EN_MTL.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.EN_MTL.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.EN_MTL.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.EN_MTL),mytitle('EN_MTL'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.EN_MTL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
% return

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
% return
%% evaluate test performance
Ypred.EN_STL = Xtest*W.EN_STL;
Ypred.EN_MTL = Xtest*W.EN_MTL;
Ypred.GN_STL = Xtest*W.GN_STL;
Ypred.GN_MTL = Xtest*W.GN_MTL;

MSE.EN_STL = norm(Ytest(:) - Ypred.EN_STL(:));
MSE.EN_MTL = norm(Ytest(:) - Ypred.EN_MTL(:));
MSE.GN_STL = norm(Ytest(:) - Ypred.GN_STL(:));
MSE.GN_MTL = norm(Ytest(:) - Ypred.GN_MTL(:));

CORR.EN_STL=corr(Ytest(:), Ypred.EN_STL(:));
CORR.EN_MTL=corr(Ytest(:), Ypred.EN_MTL(:));
CORR.GN_STL=corr(Ytest(:), Ypred.GN_STL(:));
CORR.GN_MTL=corr(Ytest(:), Ypred.GN_MTL(:));

MSE
CORR
% 
% drawnow