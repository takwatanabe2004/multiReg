%% july06_6dMNI_LADMM_FL.m
% (07/05/2014)
%=========================================================================%
% - Test out linearizezd ADMM-FL from 6dMNI structure
%=========================================================================%
%%
clear all;
purge

% rootdir = fileparts(mfilename('fullpath'));
%% load data
GRID='Grid326'; % {'Grid326','Grid1068','WashU'}

%=========================================================================%
% anomalous edge index info
%=========================================================================% 
clusterName = 'anomCluster1';
clusterPath = [get_rootdir,'/simulation/july05_groundTruthGenerator/',clusterName,'_',GRID];

[wtrue,idx_supp,C,coord] = tak_sim_weight_6d_MNI(clusterPath);
% return
figure,imexp
subplot(121),tplot(idx_supp)
subplot(122),imconnEdge(idx_supp)
drawnow

figure,imexp
subplot(121),imconn(wtrue,1)
subplot(122),imconnYeo(wtrue,GRID,1,1)
drawnow
%%
% imconnl(wtrue)
n = 100;
ntest = 100;
%=========================================================================%
% sample data matrix: iid normal or AR model
%-------------------------------------------------------------------------%
% X = randn(n,p);
% Xtest = randn(ntest,p);
%-------------------------------------------------------------------------%
T=800; rt = 0.8; rspace=[0.99,0.99,0.99];
X = tak_sample_connectome3d_subsamp_batch(n, T, coord.NSIZE, rt, rspace,coord.slex);
Xtest = tak_sample_connectome3d_subsamp_batch(ntest, T, coord.NSIZE, rt, rspace,coord.slex);
%-------------------------------------------------------------------------%
y = X*wtrue;
ytest = Xtest*wtrue;

lam=3;
gam=1;

options.rho = 1;
options.maxiter = 1000;
options.tol = 5e-4;
options.progress = 500;
options.silence = false;
options.funcval = false;
    
%=========================================================================%
% solve via admm via pcg
%=========================================================================%
[w_PCG,output]=tak_FL_regr_ADMM_PCG(X,y,lam,gam,options,C,[],wtrue);
Ypr.FLpcg = Xtest*w_PCG;
MSE.FLpcg = norm(ytest - Ypr.FLpcg);
COR.FLpcg = corr(ytest,  Ypr.FLpcg);

%=========================================================================%
% solve via linearized ADMM
%=========================================================================%
[w_LADM,output_LADM]=tak_FL_regr_LADMM(X,y,lam,gam,options,C,wtrue);
Ypr.FL_LADM = Xtest*w_LADM;
MSE.FL_LADM = norm(ytest - Ypr.FL_LADM);
COR.FL_LADM = corr(ytest,  Ypr.FL_LADM);

%=========================================================================%
% solve via linearized ADMM with one less splitting
%=========================================================================%
[w_LADM2,output_LADM2]=tak_FL_regr_LADMM_ver2(X,y,lam,gam,options,C,wtrue);
Ypr.FL_LADM2 = Xtest*w_LADM2;
MSE.FL_LADM2 = norm(ytest - Ypr.FL_LADM2);
COR.FL_LADM2 = corr(ytest,  Ypr.FL_LADM2);
%%
purge
tak_sim_plot_alg_result(output)
tak_sim_plot_alg_result(output_LADM)
tak_sim_plot_alg_result(output_LADM2)

figure,imexp
subplot(131),imconnYeo(wtrue,[],1),CAXIS=caxis;
subplot(132),imconnYeo(w_PCG,[],1),caxis(CAXIS)
subplot(133),imconnYeo(w_LADM2,[],1),caxis(CAXIS)

% figure,imexp
% subplot(121),imconn(wtrue,1),CAXIS=caxis;
% subplot(122),imconn(w_LADM2,1),caxis(CAXIS)
MSE
COR