%% june29_FL2d_MTL_ADM_FFTvsPCG.m
% (06/29/2014)
%=========================================================================%
% - ADMM-PCG vs ADMM-FFT implementation of multitask fused lasso
% - ground truth has 2-d patch structure
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)
%%
n = 200;
nx=80;
ny=80;
p = nx*ny;
q=50;

% noise level
SIGMA=3;
%% create weight vector
%=========================================================================%
% create weight vector with intra-smooth-1d structure, that has shared
% sparsity pattern across outputs (row sparsity)
%=========================================================================%
Wtrue = zeros(nx,ny,q);
k=50;
[xx,yy]=ndgrid(1:nx,1:ny);
indx1=21:40; indy1=61:70;
indx2=54:66; indy2=21:33;
indx3=11:25; indy3=11:25;
for col=1:q
    % random amplitude
    amp = (5+5*rand(3,1)).*sign(randn(3,1));    
    Wtrue(indx1,indy1,col) = amp(1) + 0.5*randn(length(indx1),length(indy1));
    Wtrue(indx2,indy2,col) = amp(2) + 0.5*randn(length(indx2),length(indy2));
    Wtrue(indx3,indy3,col) = amp(3) + 0.5*randn(length(indx3),length(indy3));
end
% figure,imexpb,CAXIS=[-1,1*10];
% subplot(141),imcov(Wtrue(:,:,1))%,caxis(CAXIS)
% subplot(142),imcov(Wtrue(:,:,2))%,caxis(CAXIS)
% subplot(143),imcov(Wtrue(:,:,3)),%caxis(CAXIS)
% subplot(144),imcov(Wtrue(:,:,4)),%caxis(CAXIS)
% drawnow
% return

Wtrue=reshape(Wtrue,[p,q]);
%% generate data samples given the weight vector
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
options.tol = 1e-4;      % <- relative change in the primal variable
options.progress = 50;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = false;    % <- track function values (may slow alg.)
options.rho=3;

C=tak_diffmat_2d([nx,ny],0);

%-------------------------------------------------------------------------%
% spectral radiuses (needed for fista stepsize)
%-------------------------------------------------------------------------%
% tic
% tau_X=svds(X,1)^2
% toc
% tic
tau_X=tnormest(X)^2;
% toc
tau_C = tnormest(C'*C);
% return
%% STUFFS FOR FL-FFT-ADMM
%=========================================================================%
% circulant 4-d difference matrix (for FFT-FL)
%=========================================================================%
Ccirc = tak_diffmat([nx,ny],1);
[A,b]=tak_get_augMat_maskMat([nx,ny]);
% Whos Ccirc A b
% return

%-------------------------------------------------------------------------%
% inputs for FL-FFT-ADMM algorithm
%-------------------------------------------------------------------------%
graphInfo.C = Ccirc;
graphInfo.A = A;
graphInfo.b = b;
graphInfo.NSIZE = [nx,ny];
%% regression
lam_stl=33;
gam_stl=15;
lam_mtl=1222;
gam_mtl=gam_stl;

% for FL
lamFL_stl=15;
gamFL_stl=10;
lamFL_mtl=15;
gamFL_mtl=15;

% W.go=1;
%%
%-------------------------------------------------------------------------%
 % FL-FFT
 %-------------------------------------------------------------------------%
% FL_STL
options.MTL=false;
[W.FL_STL_FFT,out.FL_STL_FFT]=...
    tak_FL_regr_MTL_ADMM_FFT(X,Y,lamFL_stl,gamFL_stl,options,graphInfo,Wtrue);

% FL_MTL
options.MTL=true;
[W.FL_MTL_FFT,out.FL_MTL_FFT]=...
    tak_FL_regr_MTL_ADMM_FFT(X,Y,lamFL_mtl,gamFL_mtl,options,graphInfo,Wtrue);

%-------------------------------------------------------------------------%
 % FL-PCG
 %-------------------------------------------------------------------------%
% FL_STL
options.MTL=false;
[W.FL_STL,out.FL_STL]=tak_FL_regr_MTL_ADMM_PCG(X,Y,lamFL_stl,gamFL_stl,options,C,[],Wtrue);

% FL_MTL
options.MTL=true;
[W.FL_MTL,out.FL_MTL]=tak_FL_regr_MTL_ADMM_PCG(X,Y,lamFL_mtl,gamFL_mtl,options,C,[],Wtrue);

%-------------------------------------------------------------------------%
%  Elastic-net for comparison
 %-------------------------------------------------------------------------%
% EN-STL
options.MTL=false;
options.tau = 1/(tau_X+gam_stl); % stepsize for ista
[W.EN_STL,out.EN_STL]=tak_EN_regr_MTL_FISTA(X,Y,lam_stl,gam_stl,options,Wtrue);

% EN-MTL:
options.MTL=true;
options.tau = 1/(tau_X+gam_mtl);
[W.EN_MTL,out.EN_MTL]=tak_EN_regr_MTL_FISTA(X,Y,lam_mtl,gam_mtl,options,Wtrue);

%% plot results
if isfield(W, 'FL_STL_FFT')
    Ypred.FL_STL_FFT = Xtest*W.FL_STL_FFT;
    MSE.FL_STL_FFT = norm(Ytest(:) - Ypred.FL_STL_FFT(:));
    CORR.FL_STL_FFT=corr(Ytest(:), Ypred.FL_STL_FFT(:));
    figure,imexp
%     subplot(431),tplot(log10(out.FL_STL_FFT.fval(1:end))), title('log10(function value)')
    subplot(432),tplot(log10(out.FL_STL_FFT.wdist)), title('log10(||wtrue-west||)')
    subplot(433),tplot(log10(out.FL_STL_FFT.rel_changevec)), title('log10(wnew-wold)')
    subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
    subplot(435),imagesc(W.FL_STL_FFT),mytitle('FL_STL_FFT'),caxis(CAXIS),colorbar
    subplot(436),imagesc(abs(Wtrue-W.FL_STL_FFT)),title('|wtrue-west|'),caxis(CAXIS),colorbar
    
    subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
    subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
    subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

    subplot(4,3,10),imagesc(reshape(W.FL_STL_FFT(:,1),[nx,ny])),mytitle('FL_STL_FFT1'), caxis(CAXIS1);colorbar
    subplot(4,3,11),imagesc(reshape(W.FL_STL_FFT(:,2),[nx,ny])),mytitle('FL_STL_FFT2'), caxis(CAXIS2);colorbar
    subplot(4,3,12),imagesc(reshape(W.FL_STL_FFT(:,3),[nx,ny])),mytitle('FL_STL_FFT3'), caxis(CAXIS3);colorbar
    drawnow
end

if isfield(W, 'FL_MTL_FFT')
    Ypred.FL_MTL_FFT = Xtest*W.FL_MTL_FFT;
    MSE.FL_MTL_FFT = norm(Ytest(:) - Ypred.FL_MTL_FFT(:));
    CORR.FL_MTL_FFT=corr(Ytest(:), Ypred.FL_MTL_FFT(:));
    figure,imexp
%     subplot(431),tplot(log10(out.FL_MTL_FFT.fval(1:end))), title('log10(function value)')
    subplot(432),tplot(log10(out.FL_MTL_FFT.wdist)), title('log10(||wtrue-west||)')
    subplot(433),tplot(log10(out.FL_MTL_FFT.rel_changevec)), title('log10(wnew-wold)')
    subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
    subplot(435),imagesc(W.FL_MTL_FFT),mytitle('FL_MTL_FFT'),caxis(CAXIS),colorbar
    subplot(436),imagesc(abs(Wtrue-W.FL_MTL_FFT)),title('|wtrue-west|'),caxis(CAXIS),colorbar
    
    subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
    subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
    subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

    subplot(4,3,10),imagesc(reshape(W.FL_MTL_FFT(:,1),[nx,ny])),mytitle('FL_MTL_FFT1'), caxis(CAXIS1);colorbar
    subplot(4,3,11),imagesc(reshape(W.FL_MTL_FFT(:,2),[nx,ny])),mytitle('FL_MTL_FFT2'), caxis(CAXIS2);colorbar
    subplot(4,3,12),imagesc(reshape(W.FL_MTL_FFT(:,3),[nx,ny])),mytitle('FL_MTL_FFT3'), caxis(CAXIS3);colorbar
    drawnow
end

if isfield(W, 'FL_STL')
    Ypred.FL_STL = Xtest*W.FL_STL;
    MSE.FL_STL = norm(Ytest(:) - Ypred.FL_STL(:));
    CORR.FL_STL=corr(Ytest(:), Ypred.FL_STL(:));
    figure,imexp
%     subplot(431),tplot(log10(out.FL_STL.fval(1:end))), title('log10(function value)')
    subplot(432),tplot(log10(out.FL_STL.wdist)), title('log10(||wtrue-west||)')
    subplot(433),tplot(log10(out.FL_STL.rel_changevec)), title('log10(wnew-wold)')
    subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
    subplot(435),imagesc(W.FL_STL),mytitle('FL_STL'),caxis(CAXIS),colorbar
    subplot(436),imagesc(abs(Wtrue-W.FL_STL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
    
    subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
    subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
    subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

    subplot(4,3,10),imagesc(reshape(W.FL_STL(:,1),[nx,ny])),mytitle('FL_STL1'), caxis(CAXIS1);colorbar
    subplot(4,3,11),imagesc(reshape(W.FL_STL(:,2),[nx,ny])),mytitle('FL_STL2'), caxis(CAXIS2);colorbar
    subplot(4,3,12),imagesc(reshape(W.FL_STL(:,3),[nx,ny])),mytitle('FL_STL3'), caxis(CAXIS3);colorbar
    drawnow
end

if isfield(W, 'FL_MTL')
    Ypred.FL_MTL = Xtest*W.FL_MTL;
    MSE.FL_MTL = norm(Ytest(:) - Ypred.FL_MTL(:));
    CORR.FL_MTL=corr(Ytest(:), Ypred.FL_MTL(:));
    figure,imexp
%     subplot(431),tplot(log10(out.FL_MTL.fval(1:end))), title('log10(function value)')
    subplot(432),tplot(log10(out.FL_MTL.wdist)), title('log10(||wtrue-west||)')
    subplot(433),tplot(log10(out.FL_MTL.rel_changevec)), title('log10(wnew-wold)')
    subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
    subplot(435),imagesc(W.FL_MTL),mytitle('FL_MTL'),caxis(CAXIS),colorbar
    subplot(436),imagesc(abs(Wtrue-W.FL_MTL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
    
    subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
    subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
    subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

    subplot(4,3,10),imagesc(reshape(W.FL_MTL(:,1),[nx,ny])),mytitle('FL_MTL1'), caxis(CAXIS1);colorbar
    subplot(4,3,11),imagesc(reshape(W.FL_MTL(:,2),[nx,ny])),mytitle('FL_MTL2'), caxis(CAXIS2);colorbar
    subplot(4,3,12),imagesc(reshape(W.FL_MTL(:,3),[nx,ny])),mytitle('FL_MTL3'), caxis(CAXIS3);colorbar
    drawnow
end

if isfield(W, 'EN_STL')
    Ypred.EN_STL = Xtest*W.EN_STL;
    MSE.EN_STL = norm(Ytest(:) - Ypred.EN_STL(:));
    CORR.EN_STL=corr(Ytest(:), Ypred.EN_STL(:));
    figure,imexp
%     subplot(431),tplot(log10(out.EN_STL.fval(1:end))), title('log10(function value)')
    subplot(432),tplot(log10(out.EN_STL.wdist)), title('log10(||wtrue-west||)')
    subplot(433),tplot(log10(out.EN_STL.rel_changevec)), title('log10(wnew-wold)')
    subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
    subplot(435),imagesc(W.EN_STL),mytitle('EN_STL'),caxis(CAXIS),colorbar
    subplot(436),imagesc(abs(Wtrue-W.EN_STL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
    
    subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
    subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
    subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

    subplot(4,3,10),imagesc(reshape(W.EN_STL(:,1),[nx,ny])),mytitle('EN_STL1'), caxis(CAXIS1);colorbar
    subplot(4,3,11),imagesc(reshape(W.EN_STL(:,2),[nx,ny])),mytitle('EN_STL2'), caxis(CAXIS2);colorbar
    subplot(4,3,12),imagesc(reshape(W.EN_STL(:,3),[nx,ny])),mytitle('EN_STL3'), caxis(CAXIS3);colorbar
    drawnow
end

if isfield(W, 'EN_MTL')
    Ypred.EN_MTL = Xtest*W.EN_MTL;
    MSE.EN_MTL = norm(Ytest(:) - Ypred.EN_MTL(:));
    CORR.EN_MTL=corr(Ytest(:), Ypred.EN_MTL(:));
    figure,imexp
%     subplot(431),tplot(log10(out.EN_MTL.fval(1:end))), title('log10(function value)')
    subplot(432),tplot(log10(out.EN_MTL.wdist)), title('log10(||wtrue-west||)')
    subplot(433),tplot(log10(out.EN_MTL.rel_changevec)), title('log10(wnew-wold)')
    subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
    subplot(435),imagesc(W.EN_MTL),mytitle('EN_MTL'),caxis(CAXIS),colorbar
    subplot(436),imagesc(abs(Wtrue-W.EN_MTL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
    
    subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
    subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
    subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

    subplot(4,3,10),imagesc(reshape(W.EN_MTL(:,1),[nx,ny])),mytitle('EN_MTL1'), caxis(CAXIS1);colorbar
    subplot(4,3,11),imagesc(reshape(W.EN_MTL(:,2),[nx,ny])),mytitle('EN_MTL2'), caxis(CAXIS2);colorbar
    subplot(4,3,12),imagesc(reshape(W.EN_MTL(:,3),[nx,ny])),mytitle('EN_MTL3'), caxis(CAXIS3);colorbar
    drawnow
end

% return
%% evaluate test performance

MSE
CORR
% 
% drawnow