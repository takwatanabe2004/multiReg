%% june29_4d_weight_ARfeature_randomSupport.m
% (06/29/2014)
%=========================================================================%
% - same as june26_4d_weightVector.m, but here only the feature has spatio-
%   temporal structure, assigned via matrix normal: X~MN(0,ARtime,ARspace)
% - the weight vector has random support (no structure here)
% - the aim here is to see if FL & GN can be helpful when only the feature
%   has the assumed structure, not the weight vector
% 
% REMARKS
% - don't see any advantage from GN / FL....kinda expected....although
%   I should do a formal CV tuning to confirm tis finding
%-------------------------------------------------------------------------%
% - same structure as the simulation setup in our neuroimage paper, but
%   using a 11 x 11 grid structure as in our initial simulation setup
%=========================================================================%
%%
clear all;
purge

% randn('state',0)
% rand('state',0)
%% penalty values
lam_EN = 1;  
gam_EN = .5;
lam_GN = 5;
gam_GN = .25;
lam_FL = 0.6;
gam_FL = 0.05;
%% set model parameters
sig = 1;    % noise variance
n = 500;    % # training samples
ntest=100;  % # test samples

%-------------------------------------------------------------------------%
% graph information
%-------------------------------------------------------------------------%
nx = 11;
ny = 11;

d = nx*ny; % number of nodes
p = nchoosek(d,2); % number of correlations/edges
% figure,imexpl,tak_plot_sim_nodes2d(nx,ny)

%-------------------------------------------------------------------------%
% connectome from 1D x 2D ar model (comment out to use X~N(0,I))
%-------------------------------------------------------------------------%
AR.T  = 500;       % number of time points
AR.rt = 0.7;       % temporal correlation
AR.rx = 0.9;      % spatial correlation (x)
AR.ry = 0.8;      % spatial correlation (y)
AR.rspace = [AR.rx, AR.ry];
AR.NSIZE = [nx,ny,nx,ny];

%=========================================================================%
% create graph matrix for GN and FL
%=========================================================================%
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
%% assign ground truth weight vector (no structure here...
%=========================================================================%
% assign ground truth weight on the support
%=========================================================================%
idx_mask = randsample(p,round(p/15));
k = length(idx_mask)

%=========================================================================%
% create coefficients: unif([-5,-2} || [2,+5])
%=========================================================================%
wtrue=zeros(p,1);
bias=2;        
mag=3;
wtrue(idx_mask) = (mag*rand(k,1)+bias) .* sign(randn(k,1));
% wtrue=wtrue/50;
% return
%=========================================================================%
% sample data
%=========================================================================%
X     = tak_sample_connectome2d_batch(    n, AR.T, AR.NSIZE, AR.rt, AR.rspace);
Xtest = tak_sample_connectome2d_batch(ntest, AR.T, AR.NSIZE, AR.rt, AR.rspace);
Y = X*wtrue + sig*randn(n,1);
Ytest = Xtest*wtrue;
% figure,imexpl,hist(wtrue(idx_mask),100)
% imconnr(wtrue)
% figure,imexpb,stem(Ytest),zoom on,grid off
% figure,hist(Xtest(:),100)
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
[w.EN, output.EN]=tak_EN_regr_ADMM(X,Y,lam_EN,gam_EN,options,wtrue(:));
Ypred.EN = Xtest*w.EN;
MSE.EN =  norm(Ytest(:) - Ypred.EN(:));
CORR.EN = corr(Ytest(:),  Ypred.EN(:));

%=========================================================================%
%% GN
%=========================================================================%
[w.GN, output.GN]=tak_GN_regr_FISTA(X,Y,lam_GN,gam_GN,options,C,wtrue(:));
Ypred.GN = Xtest*w.GN;
MSE.GN =  norm(Ytest(:) - Ypred.GN(:));
CORR.GN = corr(Ytest(:),  Ypred.GN(:));

%=========================================================================%
%% FL
%=========================================================================%
% FL-ADMM
% [w_FL, output_FL]=tak_FL_regr_ADMM_PCG(X,Y,lam_FL,gam_FL,options,C,[],Wtrue(:));
[w.FL, output.FL]=tak_FL_regr_ADMM_FFT(X,Y,lam_FL,gam_FL,options,graphInfo,wtrue(:));
Ypred.FL = Xtest*w.FL;
MSE.FL =  norm(Ytest(:) - Ypred.FL(:));
CORR.FL = corr(Ytest(:),  Ypred.FL(:));

%%
purge
% figure,imexpb
% subplot(131),imconn(wtrue),CAXIS=caxis;
% subplot(132),imconn(w_RR),caxis(CAXIS)
% subplot(133),imconn(w_EN),caxis(CAXIS)

%-------------------------------------------------------------------------%
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(241), imconn(wtrue,1),CAXIS=caxis;colorbar off, title('TRUE')
subplot(242), imconn( w.EN,1),caxis(CAXIS);colorbar off, title('EN')
subplot(243), imconn( w.GN,1),caxis(CAXIS);colorbar off, title('GN')
subplot(244), imconn( w.FL,1),caxis(CAXIS);colorbar off, title('FL')
subplot(245), tplot(wtrue),title('TRUE')
subplot(246), tplot( w.EN),title('EN')
subplot(247), tplot( w.GN),title('GN')
subplot(248), tplot( w.FL),title('FL')

%-------------------------------------------------------------------------%
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(241),imconnEdge(wtrue)
subplot(242),imconnEdge(w.EN)
subplot(243),imconnEdge(w.GN)
subplot(244),imconnEdge(w.FL)

subplot(245),tplot(wtrue~=0),title('TRUE')
subplot(246),tplot(w.EN~=0),title('EN')
subplot(247),tplot(w.GN~=0),title('GN')
subplot(248),tplot(w.FL~=0),title('FL')
%-------------------------------------------------------------------------%
% prediction accuracy
MSE
CORR