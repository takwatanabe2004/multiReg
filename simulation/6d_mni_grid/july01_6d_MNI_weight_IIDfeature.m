%% july01_6d_MNI_weight_IIDfeature
% (07/01/2014)
%=========================================================================%
% - 
% - weight vector support from:
%          june30_assign_anomNodeClusters_and_visualize_sliceBYslice.m 
%-------------------------------------------------------------------------%
% - script largely based off june29_4d_weight_STL_and_MTL_ARfeature.m, 
%   but with X~N(0,sig I) (not ar spatio-temporal signal model)
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

% figure,imexpb
% subplot(121),tspy(C),axis normal
% subplot(122),tspy(C'*C)
% figure,imexpb
% subplot(121),tspy(Ccirc),axis normal
% subplot(122),tspy(Ccirc'*Ccirc)
% drawnow
% return

% imconnEdgel(wsupp)
% imconnEdger(wsupp)

w=zeros(size(wsupp));
idx=find(wsupp);
w(idx)=rand(size(idx));

% % purge
% figure,imconnYeo(w,GRID,1)
% imconnYeol(w,GRID,1)
% imconnYeor(w,GRID,1)
% figure,imconnYeo(w,GRID,1,1)
% imconnYeol(w,GRID,1,1)
% imconnYeor(w,GRID,1,1)
% % purge
% figure,imconnYeoEdge(w,GRID)
% figure,imconnYeoEdge(w,GRID,1)
% imconnYeoEdgel(w,GRID)
% imconnYeoEdgel(w,GRID,1)
% imconnYeoEdger(w,GRID)
% imconnYeoEdger(w,GRID,1)
% return
%% penalty values
lam_EN = 1;  
gam_EN = .5;
lam_GN = 1;
gam_GN = 20;
lam_FL = 0.5;
gam_FL = 0.1;
%% set model parameters
sig = 1;    % noise variance
n = 100;    % # training samples
ntest=100;  % # test samples
q = 20;     % # outputs

%-------------------------------------------------------------------------%
% connectome from 1D x 2D ar model (comment out to use X~N(0,I))
%-------------------------------------------------------------------------%
% AR.T  = 200;       % number of time points
% AR.rt = 0.0;       % temporal correlation
% AR.rx = 0.9;      % spatial correlation (x)
% AR.ry = 0.7;      % spatial correlation (y)
% AR.rspace = [AR.rx, AR.ry];
% AR.NSIZE = [nx,ny,nx,ny];

%=========================================================================%
% assign ground truth weight on the support
%=========================================================================%
Wtrue=zeros(p,q);
% Wtrue(idx_mask,:) = 2;
for k=1:q
    Wtrue(idx_supp,k) = 10+10*rand + 2*rand([length(idx_supp),1]);
end
Wtrue=Wtrue/50;
% figure,imexpb
% subplot(131),imconn(Wtrue(:,1)), CAXIS=caxis;
% subplot(132),imconn(Wtrue(:,2)), caxis(CAXIS)
% subplot(133),imconn(Wtrue(:,3)), caxis(CAXIS)

figure,imexpb
subplot(131),imconnYeo(Wtrue(:,1)), CAXIS=caxis;
subplot(132),imconnYeo(Wtrue(:,2)), caxis(CAXIS)
subplot(133),imconnYeo(Wtrue(:,3)), caxis(CAXIS)
return
%=========================================================================%
% sample data
%=========================================================================%
% if exist('AR','var')
%     X = tak_sample_connectome2d_batch(n, AR.T, AR.NSIZE, AR.rt, AR.rspace);
%     Xtest = tak_sample_connectome2d_batch(ntest, AR.T, AR.NSIZE, AR.rt, AR.rspace);
% else
    X = randn(n,p);
    Xtest = randn(ntest,p);
% end
Y = X*Wtrue + sig*randn(n,q);
Ytest = Xtest*Wtrue;
% figure,imexpb,stem(Ytest),zoom on,grid off
% return
%% set options
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 1e-3;      % <- relative change in the primal variable
options.progress = 50;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = false;    % <- track function values (may slow alg.)

% ADMM parameter
options.rho=2;

% Lipschitz constant for FISTA step-size
tau_X = tnormest(X)^2;
tau_C = tnormest(C'*C);
options.tau = 1/(tau_X+tau_C); % stepsize for GN-FISTA

%=========================================================================%
% STUFFS FOR FL-FFT-ADMM
%=========================================================================%
NSIZE = [nx,ny,nx,ny];

% circulant 4-d difference matrix 
Ccirc = tak_diffmat(NSIZE,1);

%-------------------------------------------------------------------------%
% lexico-indices of sampled points 
% (excludes diagonal and upper-triangular coordinates in 4d space)
%-------------------------------------------------------------------------%
idx_samp = tak_dvec(reshape(1:d^2,[d,d]));
[A,b]=tak_get_augMat_maskMat(NSIZE, idx_samp);
% Whos Ccirc A b
% return

%-------------------------------------------------------------------------%
% inputs for FL-FFT-ADMM algorithm
%-------------------------------------------------------------------------%
graphInfo.C = Ccirc;
graphInfo.A = A;
graphInfo.b = b;
graphInfo.NSIZE = NSIZE;
%% run estimation algorithm
%-------------------------------------------------------------------------%
% ridge regression
%-------------------------------------------------------------------------%
w_RR = tak_ridge_regression(X,Y,.5);
Ypred.RR = Xtest*w_RR;
MSE.RR =  norm(Ytest(:) - Ypred.RR(:));
CORR.RR = corr(Ytest(:),  Ypred.RR(:));
% purge

%=========================================================================%
%% EN
%=========================================================================%
% enet-stl
options.MTL=false;
[w.EN_STL, output.EN_STL]=tak_EN_regr_MTL_ADMM(X,Y,lam_EN,gam_EN,options,Wtrue(:));
Ypred.EN_STL = Xtest*w.EN_STL;
MSE.EN_STL =  norm(Ytest(:) - Ypred.EN_STL(:));
CORR.EN_STL = corr(Ytest(:),  Ypred.EN_STL(:));

% enet-mtl
options.MTL=true;
[w.EN_MTL, output.EN_MTL]=tak_EN_regr_MTL_ADMM(X,Y,lam_EN,gam_EN,options,Wtrue(:));
Ypred.EN_MTL = Xtest*w.EN_MTL;
MSE.EN_MTL =  norm(Ytest(:) - Ypred.EN_MTL(:));
CORR.EN_MTL = corr(Ytest(:),  Ypred.EN_MTL(:));
%=========================================================================%
%% GN
%=========================================================================%
% GN-FISTA-STL
options.MTL=false;
[w.GN_STL, output.GN_STL]=tak_GN_regr_MTL_FISTA(X,Y,lam_GN,gam_GN,options,C,Wtrue(:));
Ypred.GN_STL = Xtest*w.GN_STL;
MSE.GN_STL =  norm(Ytest(:) - Ypred.GN_STL(:));
CORR.GN_STL = corr(Ytest(:),  Ypred.GN_STL(:));

% GN-FISTA-MTL
options.MTL=true;
[w.GN_MTL, output.GN_MTL]=tak_GN_regr_MTL_FISTA(X,Y,lam_GN,gam_GN,options,C,Wtrue(:));
Ypred.GN_MTL = Xtest*w.GN_MTL;
MSE.GN_MTL =  norm(Ytest(:) - Ypred.GN_MTL(:));
CORR.GN_MTL = corr(Ytest(:),  Ypred.GN_MTL(:));
%=========================================================================%
%% FL
%=========================================================================%
% FL-ADMM-STL
options.MTL=false;
% [w_FL, output_FL]=tak_FL_regr_MTL_ADMM_PCG(X,Y,lam_FL,gam_FL,options,C,[],Wtrue(:));
[w.FL_STL, output.FL_STL]=tak_FL_regr_MTL_ADMM_FFT(X,Y,lam_FL,gam_FL,options,graphInfo,Wtrue(:));
Ypred.FL_STL = Xtest*w.FL_STL;
MSE.FL_STL =  norm(Ytest(:) - Ypred.FL_STL(:));
CORR.FL_STL = corr(Ytest(:),  Ypred.FL_STL(:));


% FL-ADMM-MTL
options.MTL=true;
[w.FL_MTL, output.FL_MTL]=tak_FL_regr_MTL_ADMM_FFT(X,Y,lam_FL,gam_FL,options,graphInfo,Wtrue(:));
Ypred.FL_MTL = Xtest*w.FL_MTL;
MSE.FL_MTL =  norm(Ytest(:) - Ypred.FL_MTL(:));
CORR.FL_MTL = corr(Ytest(:),  Ypred.FL_MTL(:));
%%
purge
% figure,imexpb
% subplot(131),imconn(wtrue),CAXIS=caxis;
% subplot(132),imconn(w_RR),caxis(CAXIS)
% subplot(133),imconn(w_EN),caxis(CAXIS)

%-------------------------------------------------------------------------%
% STL results
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(3,4,1), imconn(Wtrue(:,1),1),CAXIS=caxis;colorbar off, title('TRUE (STL)')
subplot(3,4,2), imconn( w.EN_STL(:,1),1),caxis(CAXIS);colorbar off, title('EN_STL')
subplot(3,4,3), imconn( w.GN_STL(:,1),1),caxis(CAXIS);colorbar off, title('GN_STL')
subplot(3,4,4), imconn( w.FL_STL(:,1),1),caxis(CAXIS);colorbar off, title('FL_STL')
subplot(3,4,5), imconn(Wtrue(:,2),1),CAXIS=caxis;colorbar off, title('TRUE')
subplot(3,4,6), imconn( w.EN_STL(:,2),1),caxis(CAXIS);colorbar off, %title('EN')
subplot(3,4,7), imconn( w.GN_STL(:,2),1),caxis(CAXIS);colorbar off, %title('GN')
subplot(3,4,8), imconn( w.FL_STL(:,2),1),caxis(CAXIS);colorbar off, %title('FL')
subplot(3,4,9), imconn(Wtrue(:,3),1),CAXIS=caxis;colorbar off, title('TRUE')
subplot(3,4,10),imconn( w.EN_STL(:,3),1),caxis(CAXIS);colorbar off, %title('EN')
subplot(3,4,11),imconn( w.GN_STL(:,3),1),caxis(CAXIS);colorbar off, %title('GN')
subplot(3,4,12),imconn( w.FL_STL(:,3),1),caxis(CAXIS);colorbar off, %title('FL')
%-------------------------------------------------------------------------%
% MTL results
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(3,4,1), imconn(Wtrue(:,1),1),CAXIS=caxis;colorbar off, title('TRUE (STL)')
subplot(3,4,2), imconn( w.EN_MTL(:,1),1),caxis(CAXIS);colorbar off, title('EN_MTL')
subplot(3,4,3), imconn( w.GN_MTL(:,1),1),caxis(CAXIS);colorbar off, title('GN_MTL')
subplot(3,4,4), imconn( w.FL_MTL(:,1),1),caxis(CAXIS);colorbar off, title('FL_MTL')
subplot(3,4,5), imconn(Wtrue(:,2),1),CAXIS=caxis;colorbar off, title('TRUE')
subplot(3,4,6), imconn( w.EN_MTL(:,2),1),caxis(CAXIS);colorbar off, %title('EN')
subplot(3,4,7), imconn( w.GN_MTL(:,2),1),caxis(CAXIS);colorbar off, %title('GN')
subplot(3,4,8), imconn( w.FL_MTL(:,2),1),caxis(CAXIS);colorbar off, %title('FL')
subplot(3,4,9), imconn(Wtrue(:,3),1),CAXIS=caxis;colorbar off, title('TRUE')
subplot(3,4,10),imconn( w.EN_MTL(:,3),1),caxis(CAXIS);colorbar off, %title('EN')
subplot(3,4,11),imconn( w.GN_MTL(:,3),1),caxis(CAXIS);colorbar off, %title('GN')
subplot(3,4,12),imconn( w.FL_MTL(:,3),1),caxis(CAXIS);colorbar off, %title('FL'

%-------------------------------------------------------------------------%
% STL results
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(221),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(222),imagesc(w.EN_STL),mytitle('EN_STL'),caxis(CAXIS),colorbar
subplot(223),imagesc(w.GN_STL),mytitle('GN_STL'),caxis(CAXIS),colorbar
subplot(224),imagesc(w.FL_STL),mytitle('FL_STL'),caxis(CAXIS),colorbar

%-------------------------------------------------------------------------%
% MTL results
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(221),imagesc(Wtrue),title('Wtrue'), CAXIS=caxis;colorbar
subplot(222),imagesc(w.EN_MTL),mytitle('EN_MTL'),caxis(CAXIS),colorbar
subplot(223),imagesc(w.GN_MTL),mytitle('GN_MTL'),caxis(CAXIS),colorbar
subplot(224),imagesc(w.FL_MTL),mytitle('FL_MTL'),caxis(CAXIS),colorbar

%-------------------------------------------------------------------------%
% STL results
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(221),imagesc(Wtrue~=0),title('Wtrue'), CAXIS=caxis;colorbar
subplot(222),imagesc(w.EN_STL~=0),mytitle('EN_STL'),caxis(CAXIS),colorbar
subplot(223),imagesc(w.GN_STL~=0),mytitle('GN_STL'),caxis(CAXIS),colorbar
subplot(224),imagesc(w.FL_STL~=0),mytitle('FL_STL'),caxis(CAXIS),colorbar
colormap(1-gray)
%-------------------------------------------------------------------------%
% MTL results
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(221),imagesc(Wtrue~=0),title('Wtrue'), CAXIS=caxis;colorbar
subplot(222),imagesc(w.EN_MTL~=0),mytitle('EN_MTL'),caxis(CAXIS),colorbar
subplot(223),imagesc(w.GN_MTL~=0),mytitle('GN_MTL'),caxis(CAXIS),colorbar
subplot(224),imagesc(w.FL_MTL~=0),mytitle('FL_MTL'),caxis(CAXIS),colorbar
colormap(1-gray)
%-------------------------------------------------------------------------%
% prediction accuracy
MSE
CORR