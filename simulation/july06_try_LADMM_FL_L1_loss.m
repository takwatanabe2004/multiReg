%% july06_try_LADMM_FL_L1_loss
% (07/06/2014)
%=========================================================================%
% - Try out L1 loss via LADMM
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)
[wtrue,idx_supp,supp_info,C] = tak_sim_weight_4d_tripartite;

nx=supp_info.nx;
ny=supp_info.ny;
p=size(wtrue,1);

% imconnl(wtrue)
n = 100;
ntest = 100;
%=========================================================================%
% sample data matrix: iid normal or AR model
%-------------------------------------------------------------------------%
% X = randn(n,p);
% Xtest = randn(ntest,p);
%-------------------------------------------------------------------------%
T=200; rt = 0.8; rspace=[0.8,0.8];
X=tak_sample_connectome2d_batch(n,T,[nx,ny],rt,rspace);
Xtest=tak_sample_connectome2d_batch(ntest,T,[nx,ny],rt,rspace);
%-------------------------------------------------------------------------%
y = X*wtrue;
ytest = Xtest*wtrue;

lam=3;
gam=1;

options.rho = 1;
options.maxiter = 1000;
options.tol = 5e-5;
options.progress = 500;
options.silence = false;
options.funcval = 1;
    
%=========================================================================%
% solve via L1 loss with LADMM
%=========================================================================%
[w_L1,output_L1]=tak_FL_regr_LADMM_L1_loss(X,y,lam,gam,options,C,wtrue);
% [w_L1,output_L1]=tak_FL_regr_LADMM_L1_loss(X,y,0,0,options,C,wtrue);
Ypr.FL_L1 = Xtest*w_L1;
MSE.FL_L1 = norm(ytest - Ypr.FL_L1);
COR.FL_L1 = corr(ytest,  Ypr.FL_L1);
tak_sim_plot_alg_result(output_L1)
figure,imexp
subplot(121),imconn(wtrue,1),CAXIS=caxis;
subplot(122),imconn(w_L1,1)%,caxis(CAXIS)
return

%=========================================================================%
% solve via admm via pcg
%=========================================================================%
[w,output]=tak_FL_regr_ADMM_PCG(X,y,lam,gam,options,C,[],wtrue);
Ypr.FLpcg = Xtest*w;
MSE.FLpcg = norm(ytest - Ypr.FLpcg);
COR.FLpcg = corr(ytest,  Ypr.FLpcg);

%=========================================================================%
% solve via linearized ADMM with one less splitting
%=========================================================================%
[w_LADM2,output_LADM2]=tak_FL_regr_LADMM_ver2(X,y,lam,gam,options,C,wtrue);
Ypr.FL_LADM2 = Xtest*w_LADM2;
MSE.FL_LADM2 = norm(ytest - Ypr.FL_LADM2);
COR.FL_LADM2 = corr(ytest,  Ypr.FL_LADM2);
%%
tak_sim_plot_alg_result(output)
tak_sim_plot_alg_result(output_LADM2)


figure,imexp
subplot(121),imconn(wtrue,1),CAXIS=caxis;
subplot(122),imconn(w,1),caxis(CAXIS)

figure,imexp
subplot(121),imconn(wtrue,1),CAXIS=caxis;
subplot(122),imconn(w_LADM2,1),caxis(CAXIS)
MSE
COR