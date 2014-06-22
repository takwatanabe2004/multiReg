%% june18_1d_GN_ADMM_vs_FISTA
% (06/18/2014)
%=========================================================================%
% - compare admm and fista implementation of grahpnet regression
%=========================================================================%
%%
clear all;
purge
randn('state',0)
n = 200;
p = 1000;

k=100;
wtrue = zeros(p,1);
wtrue(1:k) = 2*sin(.02*pi*(1:k));
wtrue(1+k:2*k) = 5*sin(.02*pi*(1:k));
wtrue(601:600+k) = 4*sin(.02*pi*(1:k));
wtrue(701:700+k) = 2*sin(.02*pi*(1:k));
% data = rand(p,1); data(data>0.8) = 1; data(data<=0.8)=0;
% data = binornd(1,0.2,[p,1]);
% data(data==1)

sig = 1;
X = randn(n,p);
y = X*wtrue + sig*randn(n,1);

ntest=500;
Xtest=randn(ntest,p);
ytest=Xtest*wtrue;

% tplot(y)
% tplot(X*wtrue)
% tplott(wtrue)
% return
%% set options for admm
options.rho=1;
options.maxiter = 500;   % <- maximum number of iterations
options.tol = 5e-5;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = false; % <- display termination condition
options.funcval = true;    % <- track function values (may slow alg.)
%% estimate weight vector
%-------------------------------------------------------------------------%
% ridge regression
%-------------------------------------------------------------------------%
w_RR = tak_ridge_regression(X,y,.05);

%=========================================================================%
% graphnet-admm
%=========================================================================%

lam = 12;  % L1 penalty weight
gam = 20; % fused lasso penalty weight
C=tak_diffmat_1d(p,0);
tic
[w_GN_ADM,out_GN_ADM]=tak_admm_GN_regr_pcg(X,y,lam,gam,options,C,[],wtrue);
time_admm=toc
%=========================================================================%
% graphnet-ISTA
%=========================================================================%
options_ISTA=options;
% options_ISTA.maxiter = 500;   % <- maximum number of iterations
% options_ISTA.tol = 5e-5;      % <- relative change in the primal variable
% options_ISTA.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
% options_ISTA.silence = false; % <- display termination condition
% options_ISTA.funcval = 1;    % <- track function values (may slow alg.)

tau1=svds(X,1)^2;
[tau2,~,flag]=eigs(C'*C,1);
% tau2=eigs(C'*C,1)
% tau2=svds(C'*C,1)
if flag % eigs didnt converge, so returns shit value...use normest
    tau2=tnormest(C'*C);
end
options_ISTA.tau = 1/(tau1+gam*tau2);

%-------------------------------------------------------------------------%
% mini comparisons on how to compute spectral norm of the laplacian C'C
% - compute spectral norm of C, then square (more efficient when #col>>#row)
% - compute max eigenvalue of C'C (...seems to be more stable)
%-------------------------------------------------------------------------%
% normest_C=normest(C)^2
% SVDS_C=svds(C,1)^2
% SVDS_CtC=svds(C'*C,1)
% EIGS=eigs(C'*C,1)
% return
tic
[w_GN_IST,out_GN_IST]=tak_GN_regre_ISTA(X,y,lam,gam,options_ISTA,C,wtrue);
time_ista=toc

%=========================================================================%
% gnet-fista
%=========================================================================%
tic
[w_GN_FIST,out_GN_FIST]=tak_GN_regre_FISTA(X,y,lam,gam,options_ISTA,C,wtrue);
time_fista=toc
%% compare
figure,imexp
subplot(241)
    tplot3(log10(out_GN_ADM.fval),log10(out_GN_IST.fval),log10(out_GN_FIST.fval))
    legend('ADMM','ISTA','FISTA'),title('log10(function value)')
subplot(242)
    tplot3(log10(out_GN_ADM.wdist),log10(out_GN_IST.wdist),log10(out_GN_FIST.wdist))
    legend('ADMM','ISTA','FISTA'),title('log10(||wtrue-west||)')
subplot(243),
    tplot3(log10(out_GN_ADM.rel_changevec),log10(out_GN_IST.rel_changevec),...
           log10(out_GN_FIST.rel_changevec)), 
    legend('ADMM','ISTA','FISTA'),title('log10(wnew-wold)')
subplot(244),tplot(wtrue)
subplot(245),tplot(w_GN_ADM)
subplot(246),tplot(w_GN_IST)
subplot(247),tplot3(w_GN_ADM,w_GN_IST,w_GN_FIST), 
    HH=legend('w_est_ADM','w_est_IST','w_est_FIST');set(HH,'Interpreter','none')
subplot(248),
    tplot(abs(w_GN_ADM-w_GN_FIST)),title('|w_ADMM-w_FISTA|','Interpreter','none')
% return
%% compare methods
% purge
figure,imexp
subplot(411),stairs(wtrue),YLIM=ylim; YLIM = [YLIM(1)-1,YLIM(2)+1];
subplot(411),tstairs2(wtrue,w_RR),ylim(YLIM)
subplot(412),tstairs2(wtrue,w_GN_ADM),ylim(YLIM)
subplot(413),tstairs2(wtrue,w_GN_IST),ylim(YLIM)
subplot(414),tstairs2(wtrue,w_GN_FIST),ylim(YLIM)
drawnow

% mse_RR = norm(ytest - Xtest*w_RR)
% mse_EN = norm(ytest - Xtest*w_EN)
% mse_GN_ADM = norm(ytest - Xtest*w_GN_ADM)
% mse_GN_IST = norm(ytest - Xtest*w_GN_IST)

% corr_RR = corr(ytest , Xtest*w_RR)
corr_GN_ADM = corr(ytest , Xtest*w_GN_ADM)
corr_GN_IST = corr(ytest , Xtest*w_GN_IST)
corr_GN_FIST = corr(ytest , Xtest*w_GN_FIST)