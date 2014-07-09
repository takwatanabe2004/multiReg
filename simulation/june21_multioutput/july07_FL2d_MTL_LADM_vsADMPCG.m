%% july07_FL2d_MTL_LADM_vsADMPCG.m
% (07/07/2014)
%=========================================================================%
% - ADMM-PCG vs LADMM implementation of multitask fused lasso
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
% FL-LADMM
%-------------------------------------------------------------------------%
% FL_STL
options.MTL=false;
[W.FL_STL_LADM,out.FL_STL_LADM]=tak_FL_regr_MTL_LADMM(X,Y,lamFL_stl,gamFL_stl,options,C,Wtrue);
Ypred.FL_STL_LADM = Xtest*W.FL_STL_LADM;
MSE.FL_STL_LADM = norm(Ytest(:) - Ypred.FL_STL_LADM(:));
CORR.FL_STL_LADM=corr(Ytest(:), Ypred.FL_STL_LADM(:));

% FL_MTL
options.MTL=true;
[W.FL_MTL_LADM,out.FL_MTL_LADM]=tak_FL_regr_MTL_LADMM(X,Y,lamFL_mtl,gamFL_mtl,options,C,Wtrue);
Ypred.FL_MTL_LADM = Xtest*W.FL_MTL_LADM;
MSE.FL_MTL_LADM = norm(Ytest(:) - Ypred.FL_MTL_LADM(:));
CORR.FL_MTL_LADM=corr(Ytest(:), Ypred.FL_MTL_LADM(:));

%-------------------------------------------------------------------------%
% FL-PCG
%-------------------------------------------------------------------------%
% FL_STL
options.MTL=false;
[W.FL_STL,out.FL_STL]=tak_FL_regr_MTL_ADMM_PCG(X,Y,lamFL_stl,gamFL_stl,options,C,[],Wtrue);

% FL_MTL
options.MTL=true;
[W.FL_MTL,out.FL_MTL]=tak_FL_regr_MTL_ADMM_PCG(X,Y,lamFL_mtl,gamFL_mtl,options,C,[],Wtrue);


tak_sim_plot_alg_result(out.FL_STL_LADM)
tak_sim_plot_alg_result(out.FL_MTL_LADM)
tak_sim_plot_alg_result(out.FL_STL)
tak_sim_plot_alg_result(out.FL_MTL)

pause
%% plot results
purge
%=========================================================================%
% FL-STL-LADMM
%=========================================================================%

figure,imexp
%     subplot(431),tplot(log10(out.FL_STL.fval(1:end))), title('log10(function value)')
subplot(432),tplot(log10(out.FL_STL_LADM.wdist)), title('log10(||wtrue-west||)')
subplot(433),tplot(log10(out.FL_STL_LADM.rel_changevec)), title('log10(wnew-wold)')
subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(435),imagesc(W.FL_STL_LADM),mytitle('FL_STL_LADM'),caxis(CAXIS),colorbar
subplot(436),imagesc(abs(Wtrue-W.FL_STL_LADM)),title('|wtrue-west|'),caxis(CAXIS),colorbar

subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

subplot(4,3,10),imagesc(reshape(W.FL_STL_LADM(:,1),[nx,ny])),mytitle('FL_STL_LADM1'), caxis(CAXIS1);colorbar
subplot(4,3,11),imagesc(reshape(W.FL_STL_LADM(:,2),[nx,ny])),mytitle('FL_STL_LADM2'), caxis(CAXIS2);colorbar
subplot(4,3,12),imagesc(reshape(W.FL_STL_LADM(:,3),[nx,ny])),mytitle('FL_STL_LADM3'), caxis(CAXIS3);colorbar
drawnow

%=========================================================================%
% FL-MTL-LADMM
%=========================================================================%
figure,imexp
%     subplot(431),tplot(log10(out.FL_MTL.fval(1:end))), title('log10(function value)')
subplot(432),tplot(log10(out.FL_MTL_LADM.wdist)), title('log10(||wtrue-west||)')
subplot(433),tplot(log10(out.FL_MTL_LADM.rel_changevec)), title('log10(wnew-wold)')
subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(435),imagesc(W.FL_MTL_LADM),mytitle('FL_MTL_LADM'),caxis(CAXIS),colorbar
subplot(436),imagesc(abs(Wtrue-W.FL_MTL_LADM)),title('|wtrue-west|'),caxis(CAXIS),colorbar

subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

subplot(4,3,10),imagesc(reshape(W.FL_MTL_LADM(:,1),[nx,ny])),mytitle('FL_MTL_LADM1'), caxis(CAXIS1);colorbar
subplot(4,3,11),imagesc(reshape(W.FL_MTL_LADM(:,2),[nx,ny])),mytitle('FL_MTL_LADM2'), caxis(CAXIS2);colorbar
subplot(4,3,12),imagesc(reshape(W.FL_MTL_LADM(:,3),[nx,ny])),mytitle('FL_MTL_LADM3'), caxis(CAXIS3);colorbar
drawnow

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


% return
%% evaluate test performance

MSE
CORR