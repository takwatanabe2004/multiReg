%% july07_try_regression_FL_OCP1
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

randn('state',2)
rand('state',0)

rootdir = fileparts(mfilename('fullpath'));

load([rootdir,'/correlatedWeight_1d_patch1.mat'],'W','idxCluster')
%%
[p,q]=size(W);
%% make data
n=1200;
ntest=100;

sig=1;
X = tak_sample_AR1d(p,0.8,n);
Y = X*W + sig*randn(n,q);

Xtest = tak_sample_AR1d(p,0.8,ntest);
Ytest = Xtest*W;
% tplott(Y)
% tplott(Y-X*W),title('noise')

Ycorr = corr(Y);
% imcovvl(Ycorr)
Ycorr_abs = abs(Ycorr);
% figure,imexp,colormap(1-gray)
% for i=1:6
%     subplot(2,3,i),imagesc(Ycorr_abs>(i*0.1)),title(num2str(i*0.1)),axis image
% end
% drawnow

outputAdjmat = abs(corr(Y))>0.3;
F = tak_adjmat2incmat(outputAdjmat);
% Whos F
% figure,imexpb
% subplot(131),tspy(F)
% subplot(132),tspy(outputAdjmat)
% subplot(133),tspy(F'*F)

% figure,imexpb
% subplot(131),imcov(F)
% subplot(132),imcov(outputAdjmat)
% subplot(133),imcov(F'*F,1)
% return
%%
lam1=5;
gam1=15;
eta1=15;
options.rho = .5;
    
options.maxiter = 500;
options.tol = 1e-3;
options.progress = 50;
options.silence = false;
options.funcval = 1;

C = tak_diffmat_1d(p,0);
%% EN
options.MTL = false;
[West.EN_STL,output.EN_STL]=tak_EN_regr_MTL_ADMM(X,Y,lam1,gam1,options,W);
Ypr.EN_STL = Xtest*West.EN_STL;
CORR.EN_STL = corr(Ytest(:), Ypr.EN_STL(:));
% 
% options.MTL = true;
% [West.EN_MTL,output.EN_MTL]=tak_EN_regr_MTL_ADMM(X,Y,lam1,gam1,options,W);
% Ypr.EN_MTL = Xtest*West.EN_MTL;
% CORR.EN_MTL = corr(Ytest(:), Ypr.EN_MTL(:));
%% FL_OCP1         [W,output]=tak_FL_OCP1_regr_LADMM(X,Y,lam,gam,eta,options,C,F,wtrue)
[West.FL_OCP1,output.FL_OCP1]=tak_FL_OCP1_regr_LADMM(X,Y,lam1,gam1,eta1,options,C,F,W);
Ypr.FL_OCP1 = Xtest*West.FL_OCP1;
CORR.FL_OCP1 = corr(Ytest(:), Ypr.FL_OCP1(:));
% return
%% FL-LADM
options.MTL = false;
[West.FL_STL,output.FL_STL]=tak_FL_regr_MTL_LADMM(X,Y,lam1,gam1,options,C,W);
Ypr.FL_STL = Xtest*West.FL_STL;
CORR.FL_STL = corr(Ytest(:), Ypr.FL_STL(:));

% options.MTL = true;
% [West.FL_MTL,output.FL_MTL]=tak_FL_regr_MTL_LADMM(X,Y,lam1,gam1,options,C,W);
% Ypr.FL_MTL = Xtest*West.FL_MTL;
% CORR.FL_MTL = corr(Ytest(:), Ypr.FL_MTL(:));
%%
purge
tak_sim_plot_alg_result(output.EN_STL),title('EN-STL')
% tak_sim_plot_alg_result(output.EN_MTL),title('EN-MTL')
tak_sim_plot_alg_result(output.FL_OCP1),title('FL-OCP1')
tak_sim_plot_alg_result(output.FL_STL),title('FL-STL')
% tak_sim_plot_alg_result(output.FL_MTL),title('FL-MTL')
CORR
return
%%
purge
figure,imexp
subplot(131),imagesc(W)
subplot(132),imagesc(West.EN_STL)
% subplot(133),imagesc(West.EN_MTL)

figure,imexp
subplot(141),imagesc(W)
subplot(142),imagesc(West.FL_OCP1), mytitle('FL_OCP1')
subplot(143),imagesc(West.FL_STL), mytitle('FL_STL')
% subplot(144),imagesc(West.FL_MTL), title('FL_MTL')

figure,imexp
subplot(131),imagesc(W~=0)
subplot(132),imagesc(West.EN_STL~=0)
% subplot(133),imagesc(West.EN_MTL~=0)

figure,imexp
subplot(141),imagesc(W~=0)
subplot(142),imagesc(West.FL_OCP1~=0), mytitle('FL_OCP1')
subplot(143),imagesc(West.FL_STL~=0), mytitle('FL_STL')
% subplot(144),imagesc(West.FL_MTL~=0), mytitle('FL_MTL')
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

% figure,imexp
% subplot(311),tplot(West.EN_MTL(:,idxCluster{1}))
% subplot(312),tplot(West.EN_MTL(:,idxCluster{2}))
% subplot(313),tplot(West.EN_MTL(:,idxCluster{3}))

figure,imexp
subplot(311),tplot(West.FL_OCP1(:,idxCluster{1}))
subplot(312),tplot(West.FL_OCP1(:,idxCluster{2}))
subplot(313),tplot(West.FL_OCP1(:,idxCluster{3}))

figure,imexp
subplot(311),tplot(West.FL_STL(:,idxCluster{1}))
subplot(312),tplot(West.FL_STL(:,idxCluster{2}))
subplot(313),tplot(West.FL_STL(:,idxCluster{3}))

% figure,imexp
% subplot(311),tplot(West.FL_MTL(:,idxCluster{1}))
% subplot(312),tplot(West.FL_MTL(:,idxCluster{2}))
% subplot(313),tplot(West.FL_MTL(:,idxCluster{3}))