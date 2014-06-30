%% june29_2d_weight_ARfeature_randomSupport.m
% (06/29/2014)
%=========================================================================%
% - same as june29_4d_weight_ARfeature_randomSupport.m, but here we assign
%   a 2-d AR model on the features from 2-D space.
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

randn('state',0)
rand('state',0)
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
% 2-d feature information
%-------------------------------------------------------------------------%
nx = 60;
ny = 60;
p = nx*ny; 

%=========================================================================%
% create graph matrix for GN and FL
%=========================================================================%
% 2d difference matrix
C = tak_diffmat([nx,ny],0);
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
idx_mask = randsample(p,round(p/8));
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
% spatial correlation in the 2d features
rx = 0.9; 
ry = 0.9; 
[~,AA]=tak_prec_ar2d(nx,ny,rx,ry);
X=(AA\randn(p,n))'; 
Xtest=(AA\randn(p,ntest))'; 
Y = X*wtrue + sig*randn(n,1);
Ytest = Xtest*wtrue;

% imcovvl(reshape(X(1,:),[nx,ny]))
% figure,imexpl,hist(wtrue(idx_mask),100)
% imcovvr(reshape(wtrue,[nx,ny]))
% imedgel(reshape(wtrue,[nx,ny]))
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
NSIZE = [nx,ny];

% circulant 4-d difference matrix 
Ccirc = tak_diffmat(NSIZE,1);

%-------------------------------------------------------------------------%
% lexico-indices of sampled points 
%-------------------------------------------------------------------------%
[A,b]=tak_get_augMat_maskMat(NSIZE);
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

wtrueMat=reshape(wtrue,[nx,ny]);
w_EN_Mat=reshape(w.EN,[nx,ny]);
w_GN_Mat=reshape(w.GN,[nx,ny]);
w_FL_Mat=reshape(w.FL,[nx,ny]);
%-------------------------------------------------------------------------%
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(241), imcov(wtrueMat),CAXIS=caxis;colorbar off, title('TRUE')
subplot(242), imcov(w_EN_Mat),caxis(CAXIS);colorbar off, title('EN')
subplot(243), imcov(w_GN_Mat),caxis(CAXIS);colorbar off, title('GN')
subplot(244), imcov(w_FL_Mat),caxis(CAXIS);colorbar off, title('FL')
subplot(245), tplot(wtrue),title('TRUE')
subplot(246), tplot( w.EN),title('EN')
subplot(247), tplot( w.GN),title('GN')
subplot(248), tplot( w.FL),title('FL')

%-------------------------------------------------------------------------%
figure,set(gcf,'Units','pixels','Position', [1 51 1920 905])
subplot(241),imedge(wtrueMat)
subplot(242),imedge(w_EN_Mat)
subplot(243),imedge(w_GN_Mat)
subplot(244),imedge(w_FL_Mat)

subplot(245),tplot(wtrue~=0),title('TRUE')
subplot(246),tplot(w.EN~=0),title('EN')
subplot(247),tplot(w.GN~=0),title('GN')
subplot(248),tplot(w.FL~=0),title('FL')
%-------------------------------------------------------------------------%
% prediction accuracy
MSE
CORR