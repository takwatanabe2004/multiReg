%% july09_try_regr_ENOCP2_ADM_vs_FISTA.m
% (07/09/2014)
%=========================================================================%
% - Test scripts on dataset created from july09_corrWeight_dictMethod.m
% - Compare EN-OCP implemented via ADMM and FISTA.
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)

rootdir = fileparts(mfilename('fullpath'));

load([rootdir,'/correlatedWeight_scattered1.mat'],'W','idxCluster','dict')

[p,q]=size(W);
%% make data
n=1200;
ntest=100;

sig=0;
X=randn(n,p);
Xtest=randn(ntest,p);
% X = tak_sample_AR1d(p,0.8,n);
% Xtest = tak_sample_AR1d(p,0.8,ntest);

Y = X*W + sig*randn(n,q);
Ytest = Xtest*W;

%-------------------------------------------------------------------------%
% see if noise is at a "reasonable" level by inspection
%-------------------------------------------------------------------------%
% figure,imexpb
% subplot(131),tplot(Y),YLIM=ylim;
% subplot(132),tplot(X*W),title('signal')
% subplot(133),tplot(Y-X*W),title('noise'),ylim(YLIM)
% return
%=========================================================================%
% create output correlation matrix
%=========================================================================%
% simple correlation
Ycorr = corr(Y);
Ycorr_abs = abs(Ycorr);

% partial correlation (comment line below to use simple correlation)
Ypcorr = tak_icov2pcorr(inv(Ycorr));
Ypcorr_abs = abs(Ypcorr);

%-------------------------------------------------------------------------%
% visualize thresholded correlation graphs
%-------------------------------------------------------------------------%
% % correlation graph
% imcovvl(Ycorr)
% figure,imexp,colormap(1-gray)
% for i=1:6
%     subplot(2,3,i),imagesc(Ycorr_abs>(i*0.1)),title(num2str(i*0.1)),axis image
% end
% figure,imexp
% for i=1:6
%     subplot(2,3,i),imcov((Ycorr_abs>(i*0.1)).*Ycorr),title(num2str(i*0.1)),axis image,colorbar off
% end
% drawnow
% 
% % partial correlation graph
% imcovvl(Ypcorr)
% figure,imexp,colormap(1-gray)
% for i=1:6
%     subplot(2,3,i),imagesc(Ypcorr_abs>(i*0.1)),title(num2str(i*0.1)),axis image
% end
% figure,imexp
% for i=1:6
%     subplot(2,3,i),imcov((Ypcorr_abs>(i*0.1)).*Ypcorr),title(num2str(i*0.1)),axis image,colorbar off
% end
% drawnow
% return
%% create graph matrix F
%-------------------------------------------------------------------------%
% use binary adjacency matrix
%-------------------------------------------------------------------------%
% outputAdjmat = abs(corr(Y))>0.2;
% outputAdjmat = outputAdjmat - logical(eye(q)); % subtract diagonal (no self-edges)
% F = tak_adjmat2incmat(outputAdjmat);

%-------------------------------------------------------------------------%
% use weighted adjacency matrix
%-------------------------------------------------------------------------%
outputAdjmat = (Ycorr_abs>(0.3)).*Ycorr;
outputAdjmat = outputAdjmat - diag(diag(outputAdjmat)); % zero off the diagonal (for adj2inc function)
F = tak_adjmat2incmat_july2014(outputAdjmat,0); % F=-F;

%-------------------------------------------------------------------------%
% visualize
%-------------------------------------------------------------------------%
% figure,imexpb
% subplot(131),imagesc(F),colorbar
% subplot(132),imcov(outputAdjmat)
% subplot(133),imcov(F'*F,1), title('F''F')
% return
%% set algorithm options


% augmented lagrangian parameter
options.rho = 2;
    
options.maxiter = 2000;
options.tol = 1e-4;
options.progress = 100;
options.silence = false;
options.funcval = 1;
%% EN-OCP2 - FISTA
lam1=5;
gam1=1;
eta1=5;
[West.EN_OCP2_FISTA,output.EN_OCP2_FISTA]=tak_EN_OCP2_regr_FISTA(X,Y,lam1,gam1,eta1,options,F,W);
Ypr.EN_OCP2_FISTA = Xtest*West.EN_OCP2_FISTA;
CORR.EN_OCP2_FISTA = corr(Ytest(:), Ypr.EN_OCP2_FISTA(:))
%% EN-OCP2 - ADMM
[West.EN_OCP2,output.EN_OCP2]=tak_EN_OCP2_regr_ADMM(X,Y,lam1,gam1,eta1,options,F,W);
Ypr.EN_OCP2 = Xtest*West.EN_OCP2;
CORR.EN_OCP2 = corr(Ytest(:), Ypr.EN_OCP2(:))
%%
purge
CORR
tak_sim_plot_alg_result(output.EN_OCP2_FISTA)
tak_sim_plot_alg_result(output.EN_OCP2)

figure,imexp
for i=1:6
    err = abs(Ypr.EN_OCP2_FISTA(:,i)-Ypr.EN_OCP2(:,i));
    subplot(2,3,i),tstem(err),title('Error plots')
end
    
figure,imexp
for i=1:6
    subplot(2,3,i),tstem(Ypr.EN_OCP2_FISTA(:,i)),title('FISTA')
end    
figure,imexp
for i=1:6
    subplot(2,3,i),tstem(Ypr.EN_OCP2(:,i)),title('ADMM')
end