%% july01_6d_MNI_weight_ARfeature_tryICHOLpcg.m
% (07/01/2014)
%=========================================================================%
% CONCUSION: using ICHOL as a preconditioner can save about 25% computation
% time...also, when solving PCG over q output variables, it it turns out to
% be slightly faster by applying a looper over each variables, rather than
% applying a single run of PCG with kron-prod'ed PCG matrix.
%=========================================================================%
%%
clear
purge

GRID='Grid326'; % {'Grid326','Grid1068','WashU'}

rootdir = fileparts(mfilename('fullpath'));

randn('state',0)
rand('state',0)
%% load ground-truth weight support and graphInfo of node parcellation
dataPath = [rootdir,'/anomNodeCluster2','_',GRID,'.mat'];
dataVars = {'wsupp', 'idx_anomNodeClusters'};
load(dataPath,dataVars{:})
idx_supp = find(wsupp);

% load graphInfo
dataPath = [get_rootdir,'/data_local/graphinfo/graph_info_',GRID,'.mat'];
dataVars = {'C','coord'};
load(dataPath,dataVars{:})

% feature info
d = prod(coord.nsamp);
p = nchoosek(d,2);

% circmat for fft-based FL
NSIZE = [coord.NSIZE, coord.NSIZE];
Ccirc = tak_diffmat(NSIZE,1);
[A,b] = tak_get_augMat_maskMat(NSIZE,coord.slex);
%% penalty values
lam_FL = 0.1*3;
gam_FL = 0.1*3;
%% set model parameters
n = 150;    % # training samples
ntest=100;  % # test samples
q = 5;     % # outputs

%-------------------------------------------------------------------------%
% connectome from 1D x 3D ar model (comment out to use X~N(0,I))
%-------------------------------------------------------------------------%
AR.T  = 200;       % number of time points
AR.rt = 0.75;       % temporal correlation
AR.rx = 0.6;      % spatial correlation (x)
AR.ry = 0.6;      % spatial correlation (y)
AR.rz = 0.6;      % spatial correlation (z)
AR.rspace = [AR.rx, AR.ry, AR.rz];
AR.NSIZE = [coord.nx,coord.ny,coord.nz];

%=========================================================================%
% assign ground truth weight on the support
%=========================================================================%
Wtrue=zeros(p,q);
for k=1:q
    Wtrue(idx_supp,k) = 3 + 3*rand + 2*rand([length(idx_supp),1]);
end
Wtrue=Wtrue/5;
% figure,imexpb
% subplot(131),imconn(Wtrue(:,1)), CAXIS=caxis;
% subplot(132),imconn(Wtrue(:,2)), caxis(CAXIS)
% subplot(133),imconn(Wtrue(:,3)), caxis(CAXIS)

% figure,imexpb
% subplot(131),imconnYeo(Wtrue(:,1)), CAXIS=caxis;
% subplot(132),imconnYeo(Wtrue(:,2)), caxis(CAXIS)
% subplot(133),imconnYeo(Wtrue(:,3)), caxis(CAXIS)
% % figure,imexpb
% % subplot(131),imconnYeo(Wtrue(:,1)); 
% %     caxis([min(Wtrue(idx_supp,1)),max(Wtrue(idx_supp,1))])
% % subplot(132),imconnYeo(Wtrue(:,2)), 
% %     caxis([min(Wtrue(idx_supp,2)),max(Wtrue(idx_supp,2))])
% % subplot(133),imconnYeo(Wtrue(:,3)), 
% %     caxis([min(Wtrue(idx_supp,3)),max(Wtrue(idx_supp,3))])
% drawnow
% zoom on
% return

%=========================================================================%
% sample data
%=========================================================================%
if exist('AR','var')
    X = tak_sample_connectome3d_subsamp_batch(n, AR.T, AR.NSIZE, AR.rt, AR.rspace,coord.slex);
    Xtest = tak_sample_connectome3d_subsamp_batch(ntest, AR.T, AR.NSIZE, AR.rt, AR.rspace,coord.slex);
else
    X = randn(n,p);
    Xtest = randn(ntest,p);
end
%% add noise
%-------------------------------------------------------------------------%
% check if noise level looks reasonable
%-------------------------------------------------------------------------%
sig = 40;    % noise variance
Y = X*Wtrue + sig*randn(n,q);
Ytest = Xtest*Wtrue;

%-------------------------------------------------------------------------%
% check if noise level looks reasonable
%-------------------------------------------------------------------------%
% figure,imexpb,stem(Ytest),zoom on,grid off,YLIM=ylim;
% testNoise=sig*randn(ntest,q);
% figure,imexpb,stem(testNoise),zoom on,grid off,ylim(YLIM)
% return
%% set options
options.maxiter = 20;   % <- maximum number of iterations
options.tol = 5e-5;      % <- relative change in the primal variable
options.progress = 5;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = false;    % <- track function values (may slow alg.)

% ADMM parameter
options.rho=2;

%-------------------------------------------------------------------------%
% inputs for FL-FFT-ADMM algorithm
%-------------------------------------------------------------------------%
graphInfo.C = Ccirc;
graphInfo.A = A;
graphInfo.b = b;
graphInfo.NSIZE = NSIZE;
%% run estimation algorithm
%% FL
% FL-ADMM-MTL
options.MTL=true;
PCG.tol = 1e-6;
PCG.maxiter = 500;
disp('--- FL-MTL CG ---')
% [w.FL_MTL, output.FL_MTL]=tak_FL_regr_MTL_ADMM_FFT(X,Y,lam_FL,gam_FL,options,graphInfo,Wtrue);
[w.FL_MTL, output.FL_MTL]=tak_FL_regr_MTL_ADMM_PCG(X,Y,lam_FL,gam_FL,options,C,PCG,Wtrue);
Ypred.FL_MTL = Xtest*w.FL_MTL;
MSE.FL_MTL =  norm(Ytest(:) - Ypred.FL_MTL(:));
CORR.FL_MTL = corr(Ytest(:),  Ypred.FL_MTL(:));

% FL-ADMM-MTL...using preconditioner
PCG.precond = ichol(C'*C+speye(p));
% return
disp('--- FL-MTL PCG ---')
% [w.FL_MTL, output.FL_MTL]=tak_FL_regr_MTL_ADMM_FFT(X,Y,lam_FL,gam_FL,options,graphInfo,Wtrue);
[w.FL_MTL_looper, output.FL_MTL_looper]=tak_FL_regr_MTL_ADMM_PCG(X,Y,lam_FL,gam_FL,options,C,PCG,Wtrue);
Ypred.FL_MTL_looper = Xtest*w.FL_MTL_looper;
MSE.FL_MTL_looper =  norm(Ytest(:) - Ypred.FL_MTL_looper(:));
CORR.FL_MTL_looper = corr(Ytest(:),  Ypred.FL_MTL_looper(:));

disp('--- FL-MTL CG looper---')
PCG.precond = [];
% [w.FL_MTL, output.FL_MTL]=tak_FL_regr_MTL_ADMM_FFT(X,Y,lam_FL,gam_FL,options,graphInfo,Wtrue);
[w.FL_MTL_CG_looper, output.FL_MTL_CG_looper]=tak_FL_regr_MTL_ADMM_PCGlooper(X,Y,lam_FL,gam_FL,options,C,PCG,Wtrue);
Ypred.FL_MTL_CG_looper = Xtest*w.FL_MTL_CG_looper;
MSE.FL_MTL_CG_looper =  norm(Ytest(:) - Ypred.FL_MTL_CG_looper(:));
CORR.FL_MTL_CG_looper = corr(Ytest(:),  Ypred.FL_MTL_CG_looper(:));

disp('--- FL-MTL PCG looper---')
PCG.precond = ichol(C'*C+speye(p));
% [w.FL_MTL, output.FL_MTL]=tak_FL_regr_MTL_ADMM_FFT(X,Y,lam_FL,gam_FL,options,graphInfo,Wtrue);
[w.FL_MTL_PCG_looper, output.FL_MTL_PCG_looper]=tak_FL_regr_MTL_ADMM_PCGlooper(X,Y,lam_FL,gam_FL,options,C,PCG,Wtrue);
Ypred.FL_MTL_PCG_looper = Xtest*w.FL_MTL_PCG_looper;
MSE.FL_MTL_PCG_looper =  norm(Ytest(:) - Ypred.FL_MTL_PCG_looper(:));
CORR.FL_MTL_PCG_looper = corr(Ytest(:),  Ypred.FL_MTL_PCG_looper(:));
%%
% lwidth=3;
% purge
% figure,imexp, 
% subplot(121),title('log10(||wtrue-west||)'), hold on, grid on
%                 plot( log10(output.FL_MTL.wdist), 'g', 'linewidth',lwidth)
%                 plot( log10(output.FL_MTL_PCG_looper.wdist), 'b', 'linewidth',lwidth)
%                 H=legend('FL-CG','FL-PCG'); set(H,'fontweight','b','fontsize',22);
% subplot(122),title('log10(wnew-wold)'), hold on, grid on
%                 plot( log10(output.FL_MTL.rel_changevec), 'r', 'linewidth',lwidth)
%                 plot( log10(output.FL_MTL_PCG_looper.rel_changevec), 'g', 'linewidth',lwidth)
%                 H=legend('FL-CG','FL-PCG'); set(H,'fontweight','b','fontsize',22);
% return
% figure
% subplot(231),tplot(log10(output.EN_MTL.wdist)), title('log10(||wtrue-west||)')
% subplot(234),tplot(log10(output.EN_MTL.rel_changevec)), title('log10(wnew-wold)')
% subplot(232),tplot(log10(output.GN_MTL.wdist)), title('log10(||wtrue-west||)')
% subplot(235),tplot(log10(output.GN_MTL.rel_changevec)), title('log10(wnew-wold)')
% subplot(233),tplot(log10(output.FL_MTL.wdist)), title('log10(||wtrue-west||)')
% subplot(236),tplot(log10(output.FL_MTL.rel_changevec)), title('log10(wnew-wold)')
% drawnow
% return
% tak_sim_plot_alg_result(output.EN_MTL)
% tak_sim_plot_alg_result(output.GN_MTL)
% tak_sim_plot_alg_result(output.FL_MTL)
% purge
% figure,imexpb
% subplot(131),imconn(wtrue),CAXIS=caxis;
% subplot(132),imconn(w_RR),caxis(CAXIS)
% subplot(133),imconn(w_EN),caxis(CAXIS)

%=========================================================================%
% % STL results
% imconnYeol(Wtrue(:,1),[],1)
% figure,imexpb
% subplot(121), imconnYeo( w.EN_STL(:,1),[],1);colorbar off,caxis(caxis/2), title('EN_STL')
% subplot(122), imconnYeo( w.GN_STL(:,1),[],1);colorbar off,caxis(caxis/2), title('GN_STL')
%-------------------------------------------------------------------------%
% figure,imexpb
% subplot(121), imconn( w.EN_STL(:,1),1);colorbar off,caxis(caxis/2), title('EN_STL')
% subplot(122), imconn( w.GN_STL(:,1),1);colorbar off,caxis(caxis/2), title('GN_STL')
%-------------------------------------------------------------------------%
% figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
% subplot(3,3,1), imconn(Wtrue(:,1),1),CAXIS=caxis;colorbar off, title('TRUE (STL)')
% subplot(3,3,2), imconn( w.EN_STL(:,1),1);colorbar off, caxis(caxis/2),title('EN_STL')
% subplot(3,3,3), imconn( w.GN_STL(:,1),1);colorbar off, caxis(caxis/2),title('GN_STL')
% % subplot(3,4,4), imconn( w.FL_STL(:,1),1),caxis(CAXIS);colorbar off, title('FL_STL')
% subplot(3,3,4), imconn(Wtrue(:,2),1),CAXIS=caxis;colorbar off, title('TRUE')
% subplot(3,3,5), imconn( w.EN_STL(:,2),1);colorbar off, caxis(caxis/2) %title('EN')
% subplot(3,3,6), imconn( w.GN_STL(:,2),1);colorbar off, caxis(caxis/2)%title('GN')
% % subplot(3,4,8), imconn( w.FL_STL(:,2),1),caxis(CAXIS);colorbar off, %title('FL')
% subplot(3,3,7), imconn(Wtrue(:,3),1),CAXIS=caxis;colorbar off, title('TRUE')
% subplot(3,3,8),imconn( w.EN_STL(:,3),1);colorbar off, caxis(caxis/2)%title('EN')
% subplot(3,3,9),imconn( w.GN_STL(:,3),1),colorbar off, caxis(caxis/2)%title('GN')
% % subplot(3,4,12),imconn( w.FL_STL(:,3),1),caxis(CAXIS);colorbar off, %title('FL')
%=========================================================================%
% MTL results
% figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
% imconnYeol( Wtrue(:,1),[],1);colorbar off,caxis(caxis/2), title('EN_MTL')
% imconnYeol( w.EN_MTL(:,1),[],1);colorbar off,caxis(caxis/2), title('EN_MTL')
% imconnYeol( w.GN_MTL(:,1),[],1);colorbar off,caxis(caxis/2), title('GN_MTL')
% imconnYeol( w.FL_MTL(:,1),[],1);colorbar off,caxis(caxis/2), title('FL_MTL')

% figure,imexpb
% subplot(131),imconnYeo( w.FL_MTL(:,1),[],1);colorbar off,caxis(caxis/2), title('FL_MTL')
% subplot(132),imconnYeo( w.FL_MTL_PCG(:,1),[],1);colorbar off,caxis(caxis/2), title('FL_MTL_PCG')
% subplot(133),imconnYeo( w.FL_MTL_PCG_looper(:,1),[],1);colorbar off,caxis(caxis/2), title('FL_MTL_PCG_looper')

%-------------------------------------------------------------------------%
% figure,imexpb
% subplot(121), imconn( w.EN_MTL(:,1),1);colorbar off,caxis(caxis/2), title('EN_MTL')
% subplot(122), imconn( w.GN_MTL(:,1),1);colorbar off,caxis(caxis/2), title('GN_MTL')
%-------------------------------------------------------------------------%
% subplot(331), imconn(Wtrue(:,1),1),CAXIS=caxis;colorbar off, title('TRUE (STL)')
% subplot(3,4,4), imconn( w.FL_MTL(:,1),1),caxis(CAXIS);colorbar off, title('FL_MTL')
% subplot(334), imconn(Wtrue(:,2),1),CAXIS=caxis;colorbar off, title('TRUE')
% subplot(335), imconn( w.EN_MTL(:,2),1);colorbar off, caxis(caxis/2) %title('EN')
% subplot(336), imconn( w.GN_MTL(:,2),1);colorbar off; caxis(caxis/2) %title('GN')
% subplot(3,4,8), imconn( w.FL_MTL(:,2),1),caxis(CAXIS);colorbar off, %title('FL')
% subplot(337), imconn(Wtrue(:,3),1),CAXIS=caxis;colorbar off, title('TRUE')
% subplot(338),imconn( w.EN_MTL(:,3),1),colorbar off, caxis(caxis/2)%title('EN')
% subplot(339),imconn( w.GN_MTL(:,3),1);colorbar off,caxis(caxis/2) %title('GN')
% subplot(3,4,12),imconn( w.FL_MTL(:,3),1),caxis(CAXIS);colorbar off, %title('FL'
% 
% %-------------------------------------------------------------------------%
% % STL results
% figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
% subplot(221),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(222),imagesc(w.EN_STL),mytitle('EN_STL'),caxis(CAXIS),colorbar
% subplot(223),imagesc(w.GN_STL),mytitle('GN_STL'),caxis(CAXIS),colorbar
% % subplot(224),imagesc(w.FL_STL),mytitle('FL_STL'),caxis(CAXIS),colorbar
% 
% %-------------------------------------------------------------------------%
% % MTL results
% figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
% subplot(221),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(222),imagesc(w.EN_MTL),mytitle('EN_MTL'),caxis(CAXIS),colorbar
% subplot(223),imagesc(w.GN_MTL),mytitle('GN_MTL'),caxis(CAXIS),colorbar
% % subplot(224),imagesc(w.FL_MTL),mytitle('FL_MTL'),caxis(CAXIS),colorbar

%-------------------------------------------------------------------------%
% % STL results
% figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
% subplot(221),imagesc(Wtrue~=0),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(222),imagesc(w.EN_STL~=0),mytitle('EN_STL'),caxis(CAXIS),colorbar
% subplot(223),imagesc(w.GN_STL~=0),mytitle('GN_STL'),caxis(CAXIS),colorbar
% % subplot(224),imagesc(w.FL_STL~=0),mytitle('FL_STL'),caxis(CAXIS),colorbar
% colormap(1-gray)
% %-------------------------------------------------------------------------%
% % MTL results
% figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
% subplot(221),imagesc(Wtrue~=0),title('Wtrue'), CAXIS=caxis;colorbar
% subplot(222),imagesc(w.EN_MTL~=0),mytitle('EN_MTL'),caxis(CAXIS),colorbar
% subplot(223),imagesc(w.GN_MTL~=0),mytitle('GN_MTL'),caxis(CAXIS),colorbar
% % subplot(224),imagesc(w.FL_MTL~=0),mytitle('FL_MTL'),caxis(CAXIS),colorbar
% colormap(1-gray)
%-------------------------------------------------------------------------%
% prediction accuracy
MSE
CORR
%%
% purge
% 
% % XX = struct2array(load('C:\Users\takanori\Documents\MATLAB\multiRegr\data_local\designMatrix_FC_Grid326.mat', 'X'));
% 
% imconnl(mean(XX)), 
% % caxis([0,1])
% CAXIS1=caxis;
% imconnl(XX(1,:)), 
% % caxis([0,1])
% CAXIS2=caxis;
% 
% imconnl(mean(X)),
% caxis(CAXIS1)
% % caxis([0,1])
% imconnl(X(1,:)),
% caxis(CAXIS2)
% % caxis([0,1])