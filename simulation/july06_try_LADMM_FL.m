%% mfileName
% (07/06/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge

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

options.rho = 2;
options.maxiter = 2000;
options.tol = 1e-4;
options.progress = 50;
options.silence = false;
options.funcval = false;
    
%=========================================================================%
% solve via admm via pcg
%=========================================================================%
[w,output]=tak_FL_regr_ADMM_PCG(X,y,lam,gam,options,C,[],wtrue);
tak_sim_plot_alg_result(output)
Ypr.FLpcg = Xtest*w;
MSE.FLpcg = norm(ytest - Ypr.FLpcg);
COR.FLpcg = corr(ytest,  Ypr.FLpcg);

%=========================================================================%
% solve via linearized ADMM
%=========================================================================%
[w_LADM,output_LADM]=tak_FL_regr_LADMM(X,y,lam,gam,options,C,wtrue);
tak_sim_plot_alg_result(output_LADM)
Ypr.FL_LADM = Xtest*w_LADM;
MSE.FL_LADM = norm(ytest - Ypr.FL_LADM);
COR.FL_LADM = corr(ytest,  Ypr.FL_LADM);
%%
figure,imexp
subplot(121),imconn(wtrue,1),CAXIS=caxis;
subplot(122),imconn(w,1),caxis(CAXIS)

figure,imexp
subplot(121),imconn(wtrue,1),CAXIS=caxis;
subplot(122),imconn(w_LADM,1),caxis(CAXIS)

MSE
COR