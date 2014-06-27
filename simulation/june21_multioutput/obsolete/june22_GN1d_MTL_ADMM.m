%% june22_GN1d_MTL_ADMM
% (06/22/2014)
%=========================================================================%
% - ADMM implementation of multitask GraphNet
% - ground truth has 1-d smooth structure
%-------------------------------------------------------------------------%
% - obsolete since graphnet is shitty for admm...fista much better here
%=========================================================================%
%%
clear all;
purge

error('COnclusion: stick with FISTA for GraphNet...the condition number for one of the updates depends on gamma...which sucks...')

% randn('state',0)
% rand('state',0)
%%
n = 200;
p = 700;
q=50;

% noise level
SIGMA=2;
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
%% set options for admm
options.rho=.05;
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 1e-4;      % <- relative change in the primal variable
options.progress = 25;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = true;    % <- track function values (may slow alg.)

C=tak_diffmat_1d(p,0);
%% regression
%=========================================================================%
% EN-STL: ADMM
%=========================================================================%
lam_stl=33;
gam_stl=5;
lam_mtl=33;
gam_mtl=gam_stl;
[W.EN_STL_ADM,out.EN_STL_ADM]=...
    tak_EN_regr_STL_ADMM(X,Y,lam_stl,gam_stl,options,Wtrue);

figure,imexp
subplot(231),tplot(log10(out.EN_STL_ADM.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.EN_STL_ADM.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.EN_STL_ADM.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.EN_STL_ADM),mytitle('EN_STL_ADM'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.EN_STL_ADM)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
% return

%=========================================================================%
% EN-MTL: ADMM
%=========================================================================%
[W.EN_MTL_ADM,out.EN_MTL_ADM]=...
    tak_EN_regr_MTL_ADMM(X,Y,lam_mtl,gam_mtl,options,Wtrue);

figure,imexp
subplot(231),tplot(log10(out.EN_MTL_ADM.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.EN_MTL_ADM.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.EN_MTL_ADM.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.EN_MTL_ADM),mytitle('EN_MTL_ADM'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.EN_MTL_ADM)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
% return

%=========================================================================%
% GN-STL: ADMM
%=========================================================================%
% [W.GN_STL_ADM,out.GN_STL_ADM]=...
%     tak_GN_regr_STL_ADMM_pcg(X,Y,lam_stl,gam_stl,options,C,[],Wtrue);
% figure,imexp
% subplot(231),tplot(log10(out.GN_STL_ADM.fval(2:end))), title('log10(function value)')
% subplot(232),tplot(log10(out.GN_STL_ADM.wdist)), title('log10(||wtrue-west||)')
% subplot(233),tplot(log10(out.GN_STL_ADM.rel_changevec)), title('log10(wnew-wold)')
% subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(235),imagesc(W.GN_STL_ADM),mytitle('GN_STL_ADM'),caxis(CAXIS),colorbar
% subplot(236),imagesc(abs(Wtrue-W.GN_STL_ADM)),title('|wtrue-west|'),caxis(CAXIS),colorbar
% drawnow

%=========================================================================%
% GN-MTL: ADMM
%=========================================================================%
[W.GN_MTL_ADM,out.GN_MTL_ADM]=...
    tak_GN_regr_MTL_ADMM_pcg(X,Y,lam_mtl,gam_mtl,options,C,[],Wtrue);
figure,imexp
subplot(231),tplot(log10(out.GN_MTL_ADM.fval(2:end))), title('log10(function value)')
subplot(232),tplot(log10(out.GN_MTL_ADM.wdist)), title('log10(||wtrue-west||)')
subplot(233),tplot(log10(out.GN_MTL_ADM.rel_changevec)), title('log10(wnew-wold)')
subplot(234),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(235),imagesc(W.GN_MTL_ADM),mytitle('GN_MTL_ADM'),caxis(CAXIS),colorbar
subplot(236),imagesc(abs(Wtrue-W.GN_MTL_ADM)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow
return
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