%% july09_try_regr_smooth1d_EN_GN_OCP2.m
% (07/09/2014)
%=========================================================================%
% - Test scripts on dataset created from july09_save_correlatedSmoothOutput1d.m
% - Compare EN-OCP2 and GN-OCP2, both implemented via FISTA.
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)

rootdir = fileparts(mfilename('fullpath'));

load([rootdir,'/correlatedWeight_july09_1d_smooth1.mat'],'W','idxCluster','dict')

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
outputAdjmat = (Ycorr_abs>(0.2)).*Ycorr;
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
C=tak_diffmat_1d(p,0);

% augmented lagrangian parameter
options.rho = 2;
    
options.maxiter = 2000;
options.tol = 1e-4;
options.progress = 100;
options.silence = false;
options.funcval = 0;
%% EN-OCP2 - FISTA
lam1=5;
gam1=10;
eta1=5;
[West.EN_OCP2,output.EN_OCP2]=tak_EN_OCP2_regr_FISTA(X,Y,lam1,gam1,eta1,options,F,W);
Ypr.EN_OCP2 = Xtest*West.EN_OCP2;
CORR.EN_OCP2 = corr(Ytest(:), Ypr.EN_OCP2(:))
%% GN-OCP2 - FISTA
lam1=25;
gam1=200;
eta1=1;
[West.GN_OCP2,output.GN_OCP2]=tak_GN_OCP2_regr_FISTA(X,Y,lam1,gam1,eta1,options,C,F,W);
Ypr.GN_OCP2 = Xtest*West.GN_OCP2;
CORR.GN_OCP2 = corr(Ytest(:), Ypr.GN_OCP2(:))
%%
purge
CORR
tak_sim_plot_alg_result(output.EN_OCP2)
tak_sim_plot_alg_result(output.GN_OCP2)
%%
purge
figure,imexp
for i=1:6
    err1 = Ytest(:,i)-Ypr.EN_OCP2(:,i);
    err2 = Ytest(:,i)-Ypr.GN_OCP2(:,i);
    subplot(2,3,i),tstem2(err1,err2),title('Error plots')
    legend('EN-OCP2','GN-OCP2')
end
% tstemm(Ypr.EN_STL(:,1)),imexpt
    
figure,imexp
subplot(131),imagesc(W),CAXIS=caxis;
subplot(132),imagesc(West.EN_OCP2),caxis(CAXIS),mytitle('EN_OCP2')
subplot(133),imagesc(West.GN_OCP2),caxis(CAXIS),mytitle('GN_OCP2')


figure,imexp,colormap(1-gray)
subplot(131),imsupp(W),CAXIS=caxis;
subplot(132),imsupp(West.EN_OCP2),caxis(CAXIS),mytitle('EN_OCP2')
subplot(133),imsupp(West.GN_OCP2),caxis(CAXIS),mytitle('GN_OCP2')

% figure,imexp
% subplot(131),imagesc(W),CAXIS=caxis;
% subplot(132),imagesc(abs(W-West.EN_OCP2)),caxis(CAXIS),mytitle('|W - West.EN_OCP2')
% subplot(133),imagesc(abs(W-West.GN_OCP2)),caxis(CAXIS),mytitle('|W - West.GN_OCP2')

figure,imexp
subplot(411),tplot(W(:,idxCluster{1})), mytitle('Truth')
subplot(412),tplot(W(:,idxCluster{2}))
subplot(413),tplot(W(:,idxCluster{3}))
subplot(414),tplot(W(:,idxCluster{4}))
figure,imexp
subplot(411),tplot(West.EN_OCP2(:,idxCluster{1})), mytitle('EN_OCP2')
subplot(412),tplot(West.EN_OCP2(:,idxCluster{2}))
subplot(413),tplot(West.EN_OCP2(:,idxCluster{3}))
subplot(414),tplot(West.EN_OCP2(:,idxCluster{4}))
figure,imexp
subplot(411),tplot(West.GN_OCP2(:,idxCluster{1})), mytitle('GN_OCP2')
subplot(412),tplot(West.GN_OCP2(:,idxCluster{2}))
subplot(413),tplot(West.GN_OCP2(:,idxCluster{3}))
subplot(414),tplot(West.GN_OCP2(:,idxCluster{4}))