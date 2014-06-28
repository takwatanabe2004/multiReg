%% june26_4d_weightVector
% (06/26/2014)
%=========================================================================%
% - Create 4-D structure on the "ground truth weight vector"
% - same structure as the simulation setup in our neuroimage paper, but
%   using a 11 x 11 grid structure as in our initial simulation setup
%=========================================================================%
%%
clear all;
purge

nx = 11;
ny = 11;

d = nx*ny; % number of nodes
p = nchoosek(d,2); % number of correlations/edges
%% create graph matrix for GN and FL
%-------------------------------------------------------------------------%
% lexico-indices of sampled points 
% (excludes diagonal and upper-triangular coordinates in 4d space)
%-------------------------------------------------------------------------%
idx_samp = tak_dvec(reshape(1:d^2,[d,d]));

% 4d adjacency matrix
adjmat = tak_adjmat_subsampled([nx,ny,nx,ny],idx_samp);
C = tak_adjmat2incmat(adjmat);
% L = C'*C;
% figure,imexpb
% subplot(131),tspy(adjmat)
% subplot(132),tspy(C)
% subplot(133),tspy(L)
% return
%% assign pair of node clusters
%-------------------------------------------------------------------------%
% two clusters of nodes
%-------------------------------------------------------------------------%
% first cluster
idx_nx1 = [6:10];
idx_ny1 = [8:10];

% 2nd cluster
idx_nx2 = [3:6]-1;
idx_ny2 = [3:7]-1;

%-------------------------------------------------------------------------%
% assign two cluster of nodes
%-------------------------------------------------------------------------%
% first cluster
[xx1,yy1]=ndgrid(idx_nx1,idx_ny1);
idx_anom1 = sub2ind([nx,ny],xx1(:),yy1(:))';

% 2nd cluster
[xx2,yy2]=ndgrid(idx_nx2,idx_ny2);
idx_anom2 = sub2ind([nx,ny],xx2(:),yy2(:))';

idx_anom = [idx_anom1, idx_anom2];

% figure,imexpr,tak_plot_sim_nodes2d(nx,ny,idx_anom)
% return

%% get indices of node cluster pairs
[idx_mask,maskVec,maskMat]=tak_nodes2edges_2d(nx,ny,idx_anom1,idx_anom2);
% imedgel(maskMat)
% imedger(tak_dvecinv(maskVec,0))
% 
% figure,imexpb
% subplot(131),tplot(idx_mask)
% subplot(132),tplot(maskVec)
% subplot(133),imedge(maskMat)
%% assign ground truth weight on the support
wtrue=zeros(p,1);
wtrue(idx_mask) = 2;
wtrue(idx_mask) = 2 + .2*randn(size(idx_mask));
wtrue=wtrue/15;
% imconnl(wtrue)
%% sample data
sig = 0;
n=500;
X = randn(n,p);
y = X*wtrue(:) + sig*randn(n,1);

ntest=200;
Xtest = randn(ntest,p);
ytest = Xtest*wtrue(:);
% tplott(y)
%% set options
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 5e-4;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = false;    % <- track function values (may slow alg.)

% ADMM parameter
options.rho=2;

% Lipschitz constant for FISTA step-size
tau_X = tnormest(X)^2;
tau_C = tnormest(C'*C);
options.tau = 1/(tau_X+tau_C); % stepsize for GN-FISTA
%% run estimation algorithm
%-------------------------------------------------------------------------%
% ridge regression
%-------------------------------------------------------------------------%
w_RR = tak_ridge_regression(X,y,.5);
Ypred.RR = Xtest*w_RR;
MSE.RR =  norm(ytest(:) - Ypred.RR(:));
CORR.RR = corr(ytest(:),  Ypred.RR(:));
% purge
%-------------------------------------------------------------------------%
% enet
%-------------------------------------------------------------------------%
lam1 = 1;  % L1 penalty weight
gam1 = 1; % fused lasso penalty weight

[w_EN, output_EN]=tak_EN_regr_ADMM(X,y,lam1,gam1,options,wtrue(:));
Ypred.EN = Xtest*w_EN;
MSE.EN =  norm(ytest(:) - Ypred.EN(:));
CORR.EN = corr(ytest(:),  Ypred.EN(:));
%%
%-------------------------------------------------------------------------%
% GN-FISTA
%-------------------------------------------------------------------------%
lam_GN = 1;
gam_GN = 20;
[w_GN, output_GN]=tak_GN_regr_FISTA(X,y,lam_GN,gam_GN,options,C,wtrue(:));
Ypred.GN = Xtest*w_GN;
MSE.GN =  norm(ytest(:) - Ypred.GN(:));
CORR.GN = corr(ytest(:),  Ypred.GN(:));
%%
%-------------------------------------------------------------------------%
% FL-ADMM (PCG)
%-------------------------------------------------------------------------%
lam_FL = 0.1;
gam_FL = 0.1;
[w_FL, output_FL]=tak_FL_regr_ADMM_PCG(X,y,lam_FL,gam_FL,options,C,[],wtrue(:));
% tplottl(output_FL.wdist)
Ypred.FL = Xtest*w_FL;
MSE.FL =  norm(ytest(:) - Ypred.FL(:));
CORR.FL = corr(ytest(:),  Ypred.FL(:));
%%
% figure,imexpb
% subplot(131),imconn(wtrue),CAXIS=caxis;
% subplot(132),imconn(w_RR),caxis(CAXIS)
% subplot(133),imconn(w_EN),caxis(CAXIS)

figure,imexp
subplot(231),imconn(wtrue),CAXIS=caxis;colorbar off
subplot(232),imconn(w_GN),caxis(CAXIS);colorbar off
subplot(233),imconn(w_FL),caxis(CAXIS);colorbar off

%-------------------------------------------------------------------------%
% prediction accuracy
MSE
CORR
% figure,imexp
subplot(234),tstem2(ytest,Ypred.GN)
subplot(235),tstem2(ytest,Ypred.FL)