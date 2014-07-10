%% july09_try_regr_EN_OCP1.m
% (07/09/2014)
%=========================================================================%
% - Test admm scripts on dataset created from july09_corrWeight_dictMethod.m
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
options.rho = 2;
    
options.maxiter = 500;
options.tol = 1e-3;
options.progress = 100;
options.silence = false;
options.funcval = 1;

%% EN
lam1=5;
gam1=1;
options.MTL = false;
[West.EN_STL,output.EN_STL]=tak_EN_regr_MTL_ADMM(X,Y,lam1,gam1,options,W);
Ypr.EN_STL = Xtest*West.EN_STL;
CORR.EN_STL = corr(Ytest(:), Ypr.EN_STL(:))

%% EN-OCP1
lam1=10;
gam1=1;
eta1=3;
[West.EN_OCP1,output.EN_OCP1]=tak_EN_OCP1_regr_ADMM(X,Y,lam1,gam1,eta1,options,F,W);
Ypr.EN_OCP1 = Xtest*West.EN_OCP1;
CORR.EN_OCP1 = corr(Ytest(:), Ypr.EN_OCP1(:))
%% EN-OCP2
lam1=5;
gam1=1;
eta1=0;
[West.EN_OCP2,output.EN_OCP2]=tak_EN_OCP2_regr_ADMM(X,Y,lam1,gam1,eta1,options,F,W);
Ypr.EN_OCP2 = Xtest*West.EN_OCP2;
CORR.EN_OCP2 = corr(Ytest(:), Ypr.EN_OCP2(:))
%%
CORR
tak_sim_plot_alg_result(output.EN_STL)
tak_sim_plot_alg_result(output.EN_OCP1)
tak_sim_plot_alg_result(output.EN_OCP2)
return
%%
purge
figure,imexp
for i=1:6
    err1 = Ytest(:,i)-Ypr.EN_STL(:,i);
    err2 = Ytest(:,i)-Ypr.EN_OCP1(:,i);
    err3 = Ytest(:,i)-Ypr.EN_OCP2(:,i);
    subplot(2,3,i),tstem3(err1,err2,err3),title('Error plots')
    legend('EN-STL','EN-OCP1','EN-OCP2')
end
    % tstemm(Ypr.EN_STL(:,1)),imexpt
    
figure,imexp
subplot(141),imagesc(W),CAXIS=caxis;
subplot(142),imagesc(West.EN_STL),caxis(CAXIS),mytitle('EN_STL')
subplot(143),imagesc(West.EN_OCP1),caxis(CAXIS),mytitle('EN_OCP1')
subplot(144),imagesc(West.EN_OCP2),caxis(CAXIS),mytitle('EN_OCP2')


figure,imexp,colormap(1-gray)
subplot(141),imsupp(W),CAXIS=caxis;
subplot(142),imsupp(West.EN_STL),caxis(CAXIS),mytitle('EN_STL')
subplot(143),imsupp(West.EN_OCP1),caxis(CAXIS),mytitle('EN_OCP1')
subplot(144),imsupp(West.EN_OCP2),caxis(CAXIS),mytitle('EN_OCP2')

figure,imexp
subplot(141),imagesc(W),CAXIS=caxis;
subplot(142),imagesc(abs(W-West.EN_STL)),caxis(CAXIS),mytitle('|W - West.EN_STL')
subplot(143),imagesc(abs(W-West.EN_OCP1)),caxis(CAXIS),mytitle('|W - West.EN_OCP1')
subplot(144),imagesc(abs(W-West.EN_OCP2)),caxis(CAXIS),mytitle('|W - West.EN_OCP2')