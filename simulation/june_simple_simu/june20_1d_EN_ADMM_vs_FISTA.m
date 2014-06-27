%% june20_1d_EN_ADMM_vs_FISTA
% (06/20/2014)
%=========================================================================%
% - compare admm and fista implementation of grahpnet regression...
% - compare against Beck & Teboulle's acceleration method vs ucla2's method
%  (also in Parikh & Boyd's review)
%=========================================================================%
%%
clear all;
purge
randn('state',0)
n = 100;
p = 2500;

k=100;

%-------------------------------------------------------------------------%
% randomly assign support
%-------------------------------------------------------------------------%
rand('state',0)
idx=randsample(p,k);
wtrue = zeros(p,1);
wtrue(idx) = randn(size(idx));
% tplott(wtrue),return

sig = 0;
X = randn(n,p);
y = X*wtrue + sig*randn(n,1);

ntest=200;
Xtest=randn(ntest,p);
ytest=Xtest*wtrue;

% tplot(y)
% tplot(X*wtrue)
% tplott(wtrue)
% return
%% set options for admm
options.rho=1.5;
options.maxiter = 1000;   % <- maximum number of iterations
options.tol = 1e-54;      % <- relative change in the primal variable
options.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.silence = 1; % <- display termination condition
options.funcval = true;    % <- track function values (may slow alg.)
%% estimate weight vector
%-------------------------------------------------------------------------%
% ridge regression
%-------------------------------------------------------------------------%
w_RR = tak_ridge_regression(X,y,.05);

%=========================================================================%
% en-admm
%=========================================================================%
lam = 1;  % L1 penalty weight
gam = 1; % L2
tic
[w_EN_ADM,out_EN_ADM]=tak_admm_EN_regr(X,y,lam,gam,options,wtrue);
% time_admm=toc
% return
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
options_ISTA.tau = 1/(tau1+gam);

%=========================================================================%
% gnet-ista
%=========================================================================%
tic
[w_EN_IST,out_EN_IST]=tak_EN_regr_ISTA(X,y,lam,gam,options_ISTA,wtrue);
% time_ista=toc
% return
%=========================================================================%
% gnet-fista
%=========================================================================%
tic
[w_EN_FIST,out_EN_FIST]=tak_EN_regr_FISTA(X,y,lam,gam,options_ISTA,wtrue);
% time_fista=toc
% return
%=========================================================================%
% gnet-fista: stanford version
%=========================================================================%
tic
[w_EN_FIST_ver2,out_EN_FIST_ver2]=tak_EN_regr_FISTA_stanford(X,y,lam,gam,options_ISTA,wtrue);
% time_fista_ver2=toc
%% compare
figure,imexp
subplot(241)
    tplot4(log10(out_EN_ADM.fval),  log10(out_EN_IST.fval),...
           log10(out_EN_FIST.fval), log10(out_EN_FIST_ver2.fval))
    HH=legend('ADMM','ISTA','FISTA','FIST_ver2'); set(HH,'Interpreter','none')
    title('log10(function value)');
subplot(242),
    tplot4(log10(out_EN_ADM.wdist),  log10(out_EN_IST.wdist),...
           log10(out_EN_FIST.wdist), log10(out_EN_FIST_ver2.wdist))
    HH=legend('ADMM','ISTA','FISTA','FIST_ver2'); set(HH,'Interpreter','none')
    title('log10(||wtrue-west||)')
subplot(243),
    tplot4(log10(out_EN_ADM.rel_changevec),  log10(out_EN_IST.rel_changevec),...
           log10(out_EN_FIST.rel_changevec), log10(out_EN_FIST_ver2.rel_changevec))
    HH=legend('ADMM','ISTA','FISTA','FIST_ver2'); set(HH,'Interpreter','none')
    title('log10(wnew-wold)')
subplot(244),tplot3(w_EN_ADM,w_EN_FIST_ver2,w_EN_FIST), 
    HH=legend('w_est_ADM','w_GN_FIST_ver2','w_est_FIST');set(HH,'Interpreter','none')
subplot(245),tplot(w_EN_FIST)
subplot(246),tplot(w_EN_FIST_ver2)
subplot(247),tplot(w_EN_ADM)
subplot(248),
    tplot(abs(wtrue-w_EN_FIST)),title('|wtrue-w_FISTA|','Interpreter','none')
% return
%% compare methods
% purge
YLIM=[min(wtrue)-.1*min(wtrue), max(wtrue)+.1*max(wtrue)];
figure,imexp
% subplot(421),stairs(wtrue),YLIM=ylim; YLIM = [YLIM(1)-0.1*YLIM(1),YLIM(2)+0.1*YLIM(2)];
subplot(421),stairs(wtrue),ylim(YLIM)
subplot(423),tstairs2(w_RR,wtrue),ylim(YLIM)
subplot(425),tstairs2(w_EN_ADM,wtrue),ylim(YLIM)
subplot(427),tstairs2(w_EN_FIST,wtrue),ylim(YLIM)
drawnow

lwid=2;
subplot(4,2,[2,4,6,8])
% subplot(221),imcov(wtrue),caxis(CAXIS); colorbar('location','northoutside')
stem(ytest,'linewidth',lwid),xlim([1,length(ytest)]), hold on
stem(Xtest*w_RR,'linewidth',lwid,'color','k')
stem(Xtest*w_EN_ADM,'linewidth',lwid,'color','g')
stem(Xtest*w_EN_FIST,'linewidth',lwid,'color','r')
H=legend('true','w_RR','w_EN_ADM','w_EN_FIST');set(H,'interpreter','none')
ylim([min(ytest),max(ytest)]),grid on
drawnow

mse_RR = norm(ytest - Xtest*w_RR)
% mse_EN_ADM = norm(ytest - Xtest*w_EN_ADM)
% mse_EN_IST = norm(ytest - Xtest*w_EN_IST)
mse_EN_FIST = norm(ytest - Xtest*w_EN_FIST)

corr_RR = corr(ytest , Xtest*w_RR)
% corr_EN_ADM = corr(ytest , Xtest*w_EN_ADM)
% corr_EN_IST = corr(ytest , Xtest*w_EN_IST)
corr_EN_FIST = corr(ytest , Xtest*w_EN_FIST)