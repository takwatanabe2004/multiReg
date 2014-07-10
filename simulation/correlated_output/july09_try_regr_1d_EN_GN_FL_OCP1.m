%% july09_try_regr_1d_EN_GN_FL_OCP1.m
% (07/09/2014)
%=========================================================================%
% - Test scripts on dataset created from july09_save_correlatedSmoothOutput1d.m
% - Compare EN-OCP1 (ADMM), GN-OCP1, and FL-OCP1.
%=========================================================================%
%%
clear all;
purge

randn('state',0)
rand('state',0)

rootdir = fileparts(mfilename('fullpath'));

% load([rootdir,'/correlatedWeight_july09_1d_smooth1.mat'],'W','idxCluster','dict')
load([rootdir,'/correlatedWeight_july09_1d_patch1.mat'],'W','idxCluster','dict')
[p,q]=size(W);
%% make data
n=1200;
ntest=100;

sig=2;
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
options.rho = 1;
    
options.maxiter = 2000;
options.tol = 1e-3;
options.progress = 100;
options.silence = false;
options.funcval = 0;
%% EN - ADMM (no output-correlation info encoded)
lam1=5;
gam1=10;
options.MTL = false;
[West.EN,output.EN]=tak_EN_regr_MTL_FISTA(X,Y,lam1,gam1,options,W);
Ypr.EN = Xtest*West.EN;
CORR.EN = corr(Ytest(:), Ypr.EN(:))
%% EN-OCP1 - ADMM
lam1=5;
gam1=10;
eta1=2;
[West.EN_OCP1,output.EN_OCP1]=tak_EN_OCP1_regr_ADMM(X,Y,lam1,gam1,eta1,options,F,W);
Ypr.EN_OCP1 = Xtest*West.EN_OCP1;
CORR.EN_OCP1 = corr(Ytest(:), Ypr.EN_OCP1(:))
%% GN - FISTA (no output-correlation info encoded)
lam1=25;
gam1=200;
options.MTL = false;
[West.GN,output.GN]=tak_GN_regr_MTL_FISTA(X,Y,lam1,gam1,options,C,W);
Ypr.GN = Xtest*West.GN;
CORR.GN = corr(Ytest(:), Ypr.GN(:))
%% GN-OCP1 - LADMM
lam1=15;
gam1=200;
eta1=1.5;
[West.GN_OCP1,output.GN_OCP1]=tak_GN_OCP1_regr_LADMM(X,Y,lam1,gam1,eta1,options,C,F,W);
Ypr.GN_OCP1 = Xtest*West.GN_OCP1;
CORR.GN_OCP1 = corr(Ytest(:), Ypr.GN_OCP1(:))
%% FL - LADMM (no output-correlation info encoded)
lam1=25;
gam1=50;
options.MTL = false;
[West.FL,output.FL]=tak_FL_regr_MTL_LADMM(X,Y,lam1,gam1,options,C,W);
Ypr.FL = Xtest*West.FL;
CORR.FL = corr(Ytest(:), Ypr.FL(:))
%% FL-OCP1 - LADMM
lam1=20;
gam1=80;
eta1=5;
[West.FL_OCP1,output.FL_OCP1]=tak_FL_OCP1_regr_LADMM(X,Y,lam1,gam1,eta1,options,C,F,W);
Ypr.FL_OCP1 = Xtest*West.FL_OCP1;
CORR.FL_OCP1 = corr(Ytest(:), Ypr.FL_OCP1(:))
%%
purge
CORR
tak_sim_plot_alg_result(output.EN);         subplot(131),mytitle('EN')
tak_sim_plot_alg_result(output.EN_OCP1);         subplot(131),mytitle('EN_OCP1')
tak_sim_plot_alg_result(output.GN);         subplot(131),mytitle('GN);')
tak_sim_plot_alg_result(output.GN_OCP1);         subplot(131),mytitle('GN_OCP1')
tak_sim_plot_alg_result(output.FL);         subplot(131),mytitle('FL')
tak_sim_plot_alg_result(output.FL_OCP1);         subplot(131),mytitle('FL_OCP1')

lwid=3;
figure,imexpb,
plot(log10(output.EN.wdist),'linewidth',lwid,'color','b'), hold on
plot(log10(output.EN_OCP1.wdist),'linewidth',lwid,'color','b','linestyle','--')
plot(log10(output.GN.wdist),'linewidth',lwid,'color','g'), hold on
plot(log10(output.GN_OCP1.wdist),'linewidth',lwid,'color','g','linestyle','--')
plot(log10(output.FL.wdist),'linewidth',lwid,'color','r'), hold on
plot(log10(output.FL_OCP1.wdist),'linewidth',lwid,'color','r','linestyle','--')
HH=legend('EN','EN-OCP1','GN','GN-OCP1','FL','FL-OCP1');
set(HH,'fontweight','b','fontsize',24)
grid on
axis tight
mytitle('log10(||wtrue-west||)')
return
%%
purge
figure,imexp
for i=1:6
    err1 = Ytest(:,i)-Ypr.EN(:,i);
    err2 = Ytest(:,i)-Ypr.EN_OCP1(:,i);
    subplot(2,3,i),tstem2(err1,err2),title('Error plots')
    legend('EN','EN-OCP1')
    YLIM{i}=ylim;
end
figure,imexp
for i=1:6
    err1 = Ytest(:,i)-Ypr.GN(:,i);
    err2 = Ytest(:,i)-Ypr.GN_OCP1(:,i);
    subplot(2,3,i),tstem2(err1,err2),title('Error plots')
    legend('GN','GN-OCP1')
    ylim(YLIM{i})
end
figure,imexp
for i=1:6
    err1 = Ytest(:,i)-Ypr.FL(:,i);
    err2 = Ytest(:,i)-Ypr.FL_OCP1(:,i);
    subplot(2,3,i),tstem2(err1,err2),title('Error plots')
    legend('FL','FL-OCP1')
    ylim(YLIM{i})
end

%%
purge
figure,imexp
subplot(131),imagesc(W),CAXIS=caxis;
subplot(132),imagesc(West.EN),caxis(CAXIS),mytitle('EN')
subplot(133),imagesc(West.EN_OCP1),caxis(CAXIS),mytitle('EN_OCP1')
figure,imexp
subplot(131),imagesc(W),CAXIS=caxis;
subplot(132),imagesc(West.GN),caxis(CAXIS),mytitle('GN')
subplot(133),imagesc(West.GN_OCP1),caxis(CAXIS),mytitle('GN_OCP1')
figure,imexp
subplot(131),imagesc(W),CAXIS=caxis;
subplot(132),imagesc(West.FL),caxis(CAXIS),mytitle('FL')
subplot(133),imagesc(West.FL_OCP1),caxis(CAXIS),mytitle('FL_OCP1')

figure,imexp,colormap(1-gray)
subplot(131),imsupp(W),CAXIS=caxis;
subplot(132),imsupp(West.EN),caxis(CAXIS),mytitle('EN')
subplot(133),imsupp(West.EN_OCP1),caxis(CAXIS),mytitle('EN_OCP1')
figure,imexp,colormap(1-gray)
subplot(131),imsupp(W),CAXIS=caxis;
subplot(132),imsupp(West.GN),caxis(CAXIS),mytitle('GN')
subplot(133),imsupp(West.GN_OCP1),caxis(CAXIS),mytitle('GN_OCP1')
figure,imexp,colormap(1-gray)
subplot(131),imsupp(W),CAXIS=caxis;
subplot(132),imsupp(West.FL),caxis(CAXIS),mytitle('FL')
subplot(133),imsupp(West.FL_OCP1),caxis(CAXIS),mytitle('FL_OCP1')
% figure,imexp
% subplot(131),imagesc(W),CAXIS=caxis;
% subplot(132),imagesc(abs(W-West.EN_OCP1)),caxis(CAXIS),mytitle('|W - West.EN_OCP1')
% subplot(133),imagesc(abs(W-West.GN_OCP1)),caxis(CAXIS),mytitle('|W - West.GN_OCP1')
%%
purge
figure,imexp
subplot(411),tplot(W(:,idxCluster{1})), mytitle('Truth'),YLIM1=ylim;
subplot(412),tplot(W(:,idxCluster{2})),YLIM2=ylim;
subplot(413),tplot(W(:,idxCluster{3})),YLIM3=ylim;
subplot(414),tplot(W(:,idxCluster{4})),YLIM4=ylim;

figure,imexp
subplot(411),tplot(West.EN(:,idxCluster{1})), mytitle('EN'),ylim(YLIM1)
subplot(412),tplot(West.EN(:,idxCluster{2})),ylim(YLIM2)
subplot(413),tplot(West.EN(:,idxCluster{3})),ylim(YLIM3)
subplot(414),tplot(West.EN(:,idxCluster{4})),ylim(YLIM4)
figure,imexp
subplot(411),tplot(West.EN_OCP1(:,idxCluster{1})), mytitle('EN_OCP1'),ylim(YLIM1)
subplot(412),tplot(West.EN_OCP1(:,idxCluster{2})),ylim(YLIM2)
subplot(413),tplot(West.EN_OCP1(:,idxCluster{3})),ylim(YLIM3)
subplot(414),tplot(West.EN_OCP1(:,idxCluster{4})),ylim(YLIM4)
figure,imexp
subplot(411),tplot(West.GN(:,idxCluster{1})), mytitle('GN'),ylim(YLIM1)
subplot(412),tplot(West.GN(:,idxCluster{2})),ylim(YLIM2)
subplot(413),tplot(West.GN(:,idxCluster{3})),ylim(YLIM3)
subplot(414),tplot(West.GN(:,idxCluster{4})),ylim(YLIM4)
figure,imexp
subplot(411),tplot(West.GN_OCP1(:,idxCluster{1})), mytitle('GN_OCP1'),ylim(YLIM1)
subplot(412),tplot(West.GN_OCP1(:,idxCluster{2})),ylim(YLIM2)
subplot(413),tplot(West.GN_OCP1(:,idxCluster{3})),ylim(YLIM3)
subplot(414),tplot(West.GN_OCP1(:,idxCluster{4})),ylim(YLIM4)

figure,imexp
subplot(411),tplot(West.FL(:,idxCluster{1})), mytitle('FL'),ylim(YLIM1)
subplot(412),tplot(West.FL(:,idxCluster{2})),ylim(YLIM2)
subplot(413),tplot(West.FL(:,idxCluster{3})),ylim(YLIM3)
subplot(414),tplot(West.FL(:,idxCluster{4})),ylim(YLIM4)
figure,imexp
subplot(411),tplot(West.FL_OCP1(:,idxCluster{1})), mytitle('FL_OCP1'),ylim(YLIM1)
subplot(412),tplot(West.FL_OCP1(:,idxCluster{2})),ylim(YLIM2)
subplot(413),tplot(West.FL_OCP1(:,idxCluster{3})),ylim(YLIM3)
subplot(414),tplot(West.FL_OCP1(:,idxCluster{4})),ylim(YLIM4)
