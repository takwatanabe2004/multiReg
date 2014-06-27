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

figure,imexpr,tak_plot_sim_nodes2d(nx,ny,idx_anom)
% return

%% get indices of node cluster pairs
[idx_mask,maskVec,maskMat]=tak_nodes2edges_2d(nx,ny,idx_anom1,idx_anom2);
% imedgel(maskMat)
% imedger(tak_dvecinv(maskVec,0))
% 
figure,imexpb
subplot(131),tplot(idx_mask)
subplot(132),tplot(maskVec)
subplot(133),imedge(maskMat)

%% assign ground truth weight on the support
wtrue=zeros(p,1);
wtrue(idx_mask) = 2;
wtrue(idx_mask) = 2 + .2*randn(size(idx_mask));
wtrue=wtrue/15;
imconnl(wtrue)
return
%% sample data
sig = 0;
n=1000;
X = randn(n,p);
y = X*wtrue(:) + sig*randn(n,1);

ntest=200;
Xtest = randn(ntest,p);
ytest = Xtest*wtrue(:);
tplott(y)
%% set options for admm
options.rho=2;
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 5e-4;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = true;    % <- track function values (may slow alg.)
%-------------------------------------------------------------------------%
% ridge regression
%-------------------------------------------------------------------------%
w_RR = tak_ridge_regression(X,y,.5);
%%
% purge
%-------------------------------------------------------------------------%
% enet
%-------------------------------------------------------------------------%
lam1 = 5;  % L1 penalty weight
gam1 = 1; % fused lasso penalty weight

[w_EN, output_EN]=tak_EN_regr_ADMM(X,y,lam1,gam1,options,wtrue(:));


figure,imexpb
subplot(131),imconn(wtrue),CAXIS=caxis;
subplot(132),imconn(w_RR),caxis(CAXIS)
subplot(133),imconn(w_EN),caxis(CAXIS)
