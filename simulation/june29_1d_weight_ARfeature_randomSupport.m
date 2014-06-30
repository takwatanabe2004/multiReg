%% june29_1d_weight_ARfeature_randomSupport.m
% (06/29/2014)
%=========================================================================%
% - same as june29_4d_weight_ARfeature_randomSupport.m, but here we assign
%   a 1-d AR model on the features from 1-D space.
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
% 1-d feature information
%-------------------------------------------------------------------------%
p = 2000; 

%=========================================================================%
% create graph matrix for GN and FL
%=========================================================================%
% 2d difference matrix
C = tak_diffmat(p,0);
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
rr = 0.95; 
[~,AA]=tak_prec_ar1d(p,rr);
X=(AA\randn(p,n))'; 
Xtest=(AA\randn(p,ntest))'; 
Y = X*wtrue + sig*randn(n,1);
Ytest = Xtest*wtrue;

% figure,imexpb,plot(X(1:5,:)')
% figure,imexpl,hist(wtrue(idx_mask),100)
% tplottl(wtrue)
% tplottr(wtrue~=0)
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
NSIZE = [p];

% circulant 4-d difference matrix 
Ccirc = tak_diffmat(NSIZE,1);

%-------------------------------------------------------------------------%
% lexico-indices of sampled points 
%-------------------------------------------------------------------------%
% [A,b]=tak_get_augMat_maskMat(NSIZE);
% Whos Ccirc A b
% return

%-------------------------------------------------------------------------%
% inputs for FL-FFT-ADMM algorithm
%-------------------------------------------------------------------------%
% graphInfo.C = Ccirc;
% graphInfo.A = A;
% graphInfo.b = b;
% graphInfo.NSIZE = NSIZE;
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
[w.FL, output.FL]=tak_FL_regr_ADMM_PCG(X,Y,lam_FL,gam_FL,options,C,[],wtrue(:));
% [w.FL, output.FL]=tak_FL_regr_ADMM_FFT(X,Y,lam_FL,gam_FL,options,graphInfo,wtrue(:));
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
subplot(241), tplot(wtrue),title('TRUE')
subplot(242), tplot( w.EN),title('EN')
subplot(243), tplot( w.GN),title('GN')
subplot(244), tplot( w.FL),title('FL')
subplot(245),tplot(wtrue~=0),title('TRUE')
subplot(246),tplot(w.EN~=0),title('EN')
subplot(247),tplot(w.GN~=0),title('GN')
subplot(248),tplot(w.FL~=0),title('FL')
%-------------------------------------------------------------------------%
% prediction accuracy
MSE
CORR