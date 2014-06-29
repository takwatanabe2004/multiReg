%% june29_4d_weightVector_FLpcg_vs_FLfft.m
% (06/29/2014)
%=========================================================================%
% Same as june26_4d_weightVector.m, but compare FL algorithm of PCG-ADMM 
% and FFT-ADMM with data augmentation
%-------------------------------------------------------------------------%
% - Create 4-D structure on the "ground truth weight vector"
% - same structure as the simulation setup in our neuroimage paper, but
%   using a 11 x 11 grid structure as in our initial simulation setup
%=========================================================================%
%%
clear all;
purge

nx = 11;
ny = 13;
NSIZE=[nx,ny,nx,ny];

d = nx*ny; % number of nodes
p = nchoosek(d,2); % number of correlations/edges
%% create graph matrix for GN and FL
%-------------------------------------------------------------------------%
% lexico-indices of sampled points 
% (excludes diagonal and upper-triangular coordinates in 4d space)
%-------------------------------------------------------------------------%
idx_samp = tak_dvec(reshape(1:d^2,[d,d]));

% 4d difference matrix matrix
adjmat = tak_adjmat_subsampled(NSIZE,idx_samp);
C = tak_adjmat2incmat(adjmat);
% L = C'*C;
% figure,imexpb
% subplot(131),tspy(adjmat)
% subplot(132),tspy(C)
% subplot(133),tspy(L)
% return

%=========================================================================%
% circulant 4-d difference matrix (for FFT-FL)
%=========================================================================%
Ccirc = tak_diffmat([nx,ny,nx,ny],1);
[A,b]=tak_get_augMat_maskMat(NSIZE,idx_samp);

%-------------------------------------------------------------------------%
% inputs for FL-FFT-ADMM algorithm
%-------------------------------------------------------------------------%
graphInfo.C = Ccirc;
graphInfo.A = A;
graphInfo.b = b;
graphInfo.NSIZE = NSIZE;
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
[idx_mask,maskVec,maskMat]=tak_nodes2edges_2clusters(nx,ny,idx_anom1,idx_anom2);
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
sig = 1;
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
%% run estimation algorithm
%%
%-------------------------------------------------------------------------%
% FL-ADMM (PCG)
%-------------------------------------------------------------------------%
lam_FL = 0.1;
gam_FL = 0.1;
[w.FL_PCG, output.FL_PCG]=tak_FL_regr_ADMM_PCG(X,y,lam_FL,gam_FL,options,C,[],wtrue(:));
tak_sim_plot_alg_result(output.FL_PCG)

Ypred.FL_PCG = Xtest*w.FL_PCG;
MSE.FL_PCG =  norm(ytest(:) - Ypred.FL_PCG(:));
CORR.FL_PCG = corr(ytest(:),  Ypred.FL_PCG(:));
% return
%%
%-------------------------------------------------------------------------%
% FL-ADMM (FFT)
%-------------------------------------------------------------------------%
[w.FL_FFT, output.FL_FFT,]=tak_FL_regr_ADMM_FFT(X,y,lam_FL,gam_FL,options,graphInfo,wtrue(:));
tak_sim_plot_alg_result(output.FL_FFT)
%% view results
figure,imexpb
subplot(131),imconn(wtrue),title('wtrue'), CAXIS=caxis;colorbar
subplot(132),imconn(w.FL_PCG),mytitle('west'),caxis(CAXIS),colorbar
subplot(133),imconn(abs(wtrue-w.FL_PCG)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow

figure,imexpb
subplot(131),imconn(wtrue),title('wtrue'), CAXIS=caxis;colorbar
subplot(132),imconn(w.FL_FFT),mytitle('west'),caxis(CAXIS),colorbar
subplot(133),imconn(abs(wtrue-w.FL_FFT)),title('|wtrue-west|'),caxis(CAXIS),colorbar
drawnow

Ypred.FL_FFT = Xtest*w.FL_FFT;
MSE.FL_FFT =  norm(ytest(:) - Ypred.FL_FFT(:));
CORR.FL_FFT = corr(ytest(:),  Ypred.FL_FFT(:));

% purge
% figure,imexpb
% subplot(131),imconn(wtrue),CAXIS=caxis;
% subplot(132),imconn(w_RR),caxis(CAXIS)
% subplot(133),imconn(w_EN),caxis(CAXIS)

figure,imexpb
subplot(131),imconn(w.FL_PCG),CAXIS=caxis,title('PCG')
subplot(132),imconn(w.FL_FFT),caxis(CAXIS),title('FFT')
subplot(133),imconn(abs(w.FL_FFT-w.FL_PCG))

%-------------------------------------------------------------------------%
% prediction accuracy
MSE
CORR
figure,imexpb
subplot(121),tstem2(ytest,Ypred.FL_PCG)
subplot(122),tstem(abs(Ypred.FL_FFT-Ypred.FL_PCG))