%% july07_try_regression_EN_FL_baseline
% (07/07/2014)
%=========================================================================%
% - Run EN-STL, EN-MTL, FL-STL, FL-MTL on the correlated weight vector.
% - as expected, the L1/L2-MTL here harms the performance, as it (wrongly)
%   rigidly groups together the weights across outputs.
% - this is exactly why we need OCP (output correlated penatly)
%=========================================================================%
%%
clear all;
purge

rootdir = fileparts(mfilename('fullpath'));

load([rootdir,'/correlatedWeight_1d_patch1.mat'],'W','idxCluster')
%%
[p,q]=size(W);
%% make data
n=1000;
ntest=100;

sig=0;
X = tak_sample_AR1d(p,0.8,n);
Y = X*W + sig*randn(n,q);

Xtest = tak_sample_AR1d(p,0.8,ntest);
Ytest = Xtest*W;
% tplott(Y)
% tplott(Y-X*W),title('noise')
%%
lam1=25;
gam1=25;
options.rho = 1;
    
options.maxiter = 500;
options.tol = 5e-4;
options.progress = inf;
options.silence = false;
options.funcval = false;
%% EN
options.MTL = false;
[West.EN_STL,output.EN_STL]=tak_EN_regr_MTL_ADMM(X,Y,lam1,gam1,options,W);
Ypr.EN_STL = Xtest*West.EN_STL;
CORR.EN_STL = corr(Ytest(:), Ypr.EN_STL(:));


options.MTL = true;
[West.EN_MTL,output.EN_MTL]=tak_EN_regr_MTL_ADMM(X,Y,lam1,gam1,options,W);
Ypr.EN_MTL = Xtest*West.EN_MTL;
CORR.EN_MTL = corr(Ytest(:), Ypr.EN_MTL(:));
%% FL
lam1=5;
gam1=15;
C = tak_diffmat_1d(p,0);
options.MTL = false;
[West.FL_STL,output.FL_STL]=tak_FL_regr_MTL_ADMM_MY_PCG(X,Y,lam1,gam1,options,C,[],W);
Ypr.FL_STL = Xtest*West.FL_STL;
CORR.FL_STL = corr(Ytest(:), Ypr.FL_STL(:));

options.MTL = true;
[West.FL_MTL,output.FL_MTL]=tak_FL_regr_MTL_ADMM_MY_PCG(X,Y,lam1,gam1,options,C,[],W);
Ypr.FL_MTL = Xtest*West.FL_MTL;
CORR.FL_MTL = corr(Ytest(:), Ypr.FL_MTL(:));

%%

%%
purge
tak_sim_plot_alg_result(output.EN_STL),title('EN-STL')
tak_sim_plot_alg_result(output.EN_MTL),title('EN-MTL')
tak_sim_plot_alg_result(output.FL_STL),title('FL-STL')
tak_sim_plot_alg_result(output.FL_MTL),title('FL-MTL')
CORR
return
%%
purge
figure,imexp
subplot(131),imagesc(W)
subplot(132),imagesc(West.EN_STL)
subplot(133),imagesc(West.EN_MTL)

figure,imexp
subplot(131),imagesc(W)
subplot(132),imagesc(West.FL_STL)
subplot(133),imagesc(West.FL_MTL)

figure,imexp
subplot(131),imagesc(W~=0)
subplot(132),imagesc(West.EN_STL~=0)
subplot(133),imagesc(West.EN_MTL~=0)

figure,imexp
subplot(131),imagesc(W~=0)
subplot(132),imagesc(West.FL_STL~=0)
subplot(133),imagesc(West.FL_MTL~=0)
%%
purge
figure,imexp
subplot(311),tplot(W(:,idxCluster{1}))
subplot(312),tplot(W(:,idxCluster{2}))
subplot(313),tplot(W(:,idxCluster{3}))

figure,imexp
subplot(311),tplot(West.EN_STL(:,idxCluster{1}))
subplot(312),tplot(West.EN_STL(:,idxCluster{2}))
subplot(313),tplot(West.EN_STL(:,idxCluster{3}))

figure,imexp
subplot(311),tplot(West.EN_MTL(:,idxCluster{1}))
subplot(312),tplot(West.EN_MTL(:,idxCluster{2}))
subplot(313),tplot(West.EN_MTL(:,idxCluster{3}))

figure,imexp
subplot(311),tplot(West.FL_STL(:,idxCluster{1}))
subplot(312),tplot(West.FL_STL(:,idxCluster{2}))
subplot(313),tplot(West.FL_STL(:,idxCluster{3}))

figure,imexp
subplot(311),tplot(West.FL_MTL(:,idxCluster{1}))
subplot(312),tplot(West.FL_MTL(:,idxCluster{2}))
subplot(313),tplot(West.FL_MTL(:,idxCluster{3}))