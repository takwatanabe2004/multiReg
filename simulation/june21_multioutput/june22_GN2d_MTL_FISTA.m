%% june22_GN2d_MTL_FISTA
% (06/22/2014)
%=========================================================================%
% - FISTA implementation of multitask GraphNet
% - ground truth has 2-d smooth structure
%=========================================================================%
%%
clear all;
purge

% randn('state',0)
% rand('state',0)
%%
n = 200;
nx=80;
ny=80;
p = nx*ny;
q=10;

% noise level
SIGMA=1;
%% create weight vector
%=========================================================================%
% create weight vector with intra-smooth-1d structure, that has shared
% sparsity pattern across outputs (row sparsity)
%=========================================================================%
Wtrue = zeros(nx,ny,q);
k=50;
[xx,yy]=ndgrid(1:nx,1:ny);
ctr1 = [21 61]; rad1 = 10;
ctr2 = [66 66]; rad2 = 8;
ctr3 = [44 25]; rad3 = 18;
quad_kern = @(crd) (1-crd.^2).^2 .* (abs(crd)<=1);
crd = @(ctr,rad) sqrt( (xx-ctr(1)).^2 + (yy-ctr(2)).^2 )/rad;
% crd = @(ctr,rad) sqrt( ((xx-ctr(1))/rad).^2 + ((yy-ctr(2))/rad).^2 );
crd2 = @(ctr,rad) sqrt( ((xx-ctr(1))/rad).^2 + ((yy-ctr(2))/rad).^2  + ...
                        1*(xx-ctr(1)).*(yy-ctr(2))/rad^2               );
crd3 = @(ctr,rad) sqrt( ((xx-ctr(1))/rad).^2 + ((yy-ctr(2))/rad).^2  + ...
                        -1*(xx-ctr(1)).*(yy-ctr(2))/rad^2               );
for col=1:q
    % random amplitude
    amp = (5+5*rand(3,1)).*sign(randn(3,1));    
    Wtrue(:,:,col) = ...
              + amp(1)*quad_kern(crd2(ctr1,rad1)) ...
              + amp(2)*quad_kern(crd(ctr2,rad2)) ...
              + amp(3)*quad_kern(crd3(ctr3,rad3));
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
options.tol = 1e-3;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = true;    % <- track function values (may slow alg.)

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

W.go=1;

% EN-STL
options.MTL=false;
options.tau = 1/(tau_X+gam_stl); % stepsize for ista
[W.EN_STL,out.EN_STL]=tak_EN_regr_MTL_FISTA(X,Y,lam_stl,gam_stl,options,Wtrue);

% EN-MTL:
options.MTL=true;
options.tau = 1/(tau_X+gam_mtl);
[W.EN_MTL,out.EN_MTL]=tak_EN_regr_MTL_FISTA(X,Y,lam_mtl,gam_mtl,options,Wtrue);

% GN-STL
options.MTL=false;
options.tau = 1/(tau_X+gam_stl*tau_C);
[W.GN_STL,out.GN_STL]=tak_GN_regr_MTL_FISTA(X,Y,lam_stl,5*gam_stl,options,C,Wtrue);

% GN-MTL
options.MTL=true;
options.tau = 1/(tau_X+gam_mtl*tau_C);
[W.GN_MTL,out.GN_MTL]=tak_GN_regr_MTL_FISTA(X,Y,gam_mtl,gam_mtl,options,C,Wtrue);
%% plot results
if isfield(W, 'EN_STL')
    Ypred.EN_STL = Xtest*W.EN_STL;
    MSE.EN_STL = norm(Ytest(:) - Ypred.EN_STL(:));
    CORR.EN_STL=corr(Ytest(:), Ypred.EN_STL(:));
    figure,imexp
    subplot(431),tplot(log10(out.EN_STL.fval(1:end))), title('log10(function value)')
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
    subplot(431),tplot(log10(out.EN_MTL.fval(1:end))), title('log10(function value)')
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

if isfield(W, 'GN_STL')
    Ypred.GN_STL = Xtest*W.GN_STL;
    MSE.GN_STL = norm(Ytest(:) - Ypred.GN_STL(:));
    CORR.GN_STL=corr(Ytest(:), Ypred.GN_STL(:));
    figure,imexp
    subplot(431),tplot(log10(out.GN_STL.fval(1:end))), title('log10(function value)')
    subplot(432),tplot(log10(out.GN_STL.wdist)), title('log10(||wtrue-west||)')
    subplot(433),tplot(log10(out.GN_STL.rel_changevec)), title('log10(wnew-wold)')
    subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
    subplot(435),imagesc(W.GN_STL),mytitle('GN_STL'),caxis(CAXIS),colorbar
    subplot(436),imagesc(abs(Wtrue-W.GN_STL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
    
    subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
    subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
    subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

    subplot(4,3,10),imagesc(reshape(W.GN_STL(:,1),[nx,ny])),mytitle('GN_STL1'), caxis(CAXIS1);colorbar
    subplot(4,3,11),imagesc(reshape(W.GN_STL(:,2),[nx,ny])),mytitle('GN_STL2'), caxis(CAXIS2);colorbar
    subplot(4,3,12),imagesc(reshape(W.GN_STL(:,3),[nx,ny])),mytitle('GN_STL3'), caxis(CAXIS3);colorbar
    drawnow
end

if isfield(W, 'GN_MTL')
    Ypred.GN_MTL = Xtest*W.GN_MTL;
    MSE.GN_MTL = norm(Ytest(:) - Ypred.GN_MTL(:));
    CORR.GN_MTL=corr(Ytest(:), Ypred.GN_MTL(:));
    figure,imexp
    subplot(431),tplot(log10(out.GN_MTL.fval(1:end))), title('log10(function value)')
    subplot(432),tplot(log10(out.GN_MTL.wdist)), title('log10(||wtrue-west||)')
    subplot(433),tplot(log10(out.GN_MTL.rel_changevec)), title('log10(wnew-wold)')
    subplot(434),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
    subplot(435),imagesc(W.GN_MTL),mytitle('GN_MTL'),caxis(CAXIS),colorbar
    subplot(436),imagesc(abs(Wtrue-W.GN_MTL)),title('|wtrue-west|'),caxis(CAXIS),colorbar
    
    subplot(437),imagesc(reshape(Wtrue(:,1),[nx,ny])),title('Wtrue1'), CAXIS1=caxis;colorbar
    subplot(438),imagesc(reshape(Wtrue(:,2),[nx,ny])),title('Wtrue2'), CAXIS2=caxis;colorbar
    subplot(439),imagesc(reshape(Wtrue(:,3),[nx,ny])),title('Wtrue3'), CAXIS3=caxis;colorbar

    subplot(4,3,10),imagesc(reshape(W.GN_MTL(:,1),[nx,ny])),mytitle('GN_MTL1'), caxis(CAXIS1);colorbar
    subplot(4,3,11),imagesc(reshape(W.GN_MTL(:,2),[nx,ny])),mytitle('GN_MTL2'), caxis(CAXIS2);colorbar
    subplot(4,3,12),imagesc(reshape(W.GN_MTL(:,3),[nx,ny])),mytitle('GN_MTL3'), caxis(CAXIS3);colorbar
    drawnow
end
% return
%% evaluate test performance

MSE
CORR
% 
% drawnow