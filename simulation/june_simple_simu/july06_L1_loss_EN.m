%% july06_L1_loss_EN
% (07/06/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge

p = 200;
k = round(p/10);

idx_supp = randsample(p,k);

wtrue=zeros(p,1);
wtrue(idx_supp) = tak_sample_signed_unif([4,12],k);
% tplott(wtrue)

sig=1;

n = 500;
X = randn(n,p);
y = X*wtrue + sig*randn(n,1);

% add outliers
yclean=y;
idx_outlier=randsample(n,round(n/5));
y(idx_outlier) = y(idx_outlier) + 5*tak_sample_signed_unif([4,6],length(idx_outlier));
figure,imexpt
subplot(131),tplot(yclean)
subplot(132),tplot(y)
subplot(133),tstem(yclean-y)

% return
%% ADMM
lam=1;
gam=1;

options.rho = 1;
    
options.maxiter = 500;
options.tol = 5e-4;
options.progress = inf;
options.silence = false;
options.funcval = true;

[w.L2_EN, output.L2_EN] = tak_EN_regr_ADMM(X,y,lam,gam,options,wtrue);


[w.L1_EN, output.L1_EN] = tak_EN_L1regr_ADMM(X,y,lam,gam,options,wtrue);

tak_sim_plot_alg_result(output.L2_EN)
tak_sim_plot_alg_result(output.L1_EN)
figure,imexpt,tstem3(wtrue,w.L2_EN,w.L1_EN),legend('wtrue','L2EN','L1EN')
figure,imexpt,tstem2(abs(wtrue-w.L2_EN),abs(wtrue-w.L1_EN)),legend('errL2EN','errL1EN')
figure,imexpt,tplot(w.L2_EN-w.L1_EN)