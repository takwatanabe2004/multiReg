%% june08_synthetic_data_cppd.m 
% (06/08/2014)
%=========================================================================%
% - create synthetic dataset for regression
% - use synthetic data to test my CPPD alg (Chamboll-Pock Prima-dual)
%=========================================================================%
%%
clear all
purge
load('Results_tak.mat', 'ConnMean', 'ConnVar')
load('graph_info347_2d.mat', 'coord','adjmat')
load('sim_anom_node_info_8nn.mat', 'anom_nodes')
%% set options
% set SNR devel (same as cohen's d) http://en.wikipedia.org/wiki/Cohen%27s_d#Cohen.27s_d 
% http://www.3rs-reduction.co.uk/html/6__power_and_sample_size.html
ntr = 500;
nts = 100;

p = nchoosek(coord.num_nodes,2);
%% create dataset
randn('seed',0)
%=========================================================================%
% create design matrix
%=========================================================================%
Xtr = randn(ntr,p);
Xts = randn(nts,p);

%=========================================================================%
% create ground truth weight vector
%=========================================================================%
wtrue = zeros(p,1);
wtrue(anom_nodes.idx_conn)= 0.5 + randn(size(anom_nodes.idx_conn));
% w= 10 + rand(p,1);

% normalize weight vector
% w=w/norm(w);
% norm(w)
% return

% figure,imexpb
% subplot(121),tplot(w)
% subplot(122),imcov(tak_dvecinv(w,0))
% return

%=========================================================================%
% additive noise model
%=========================================================================%
snr=.1;
ytr = Xtr*wtrue + snr*randn(ntr,1);
yts = Xts*wtrue + snr*randn(nts,1);
yts_clean = Xts*wtrue;
%% apply regression
%=========================================================================%
% set options
%=========================================================================%
% options.lambda=0;  % L1 penalty weight
% options.gamma =5; % fused lasso penalty weight

options.lambda=2^-3.5;  % L1 penalty weight
options.gamma =2^-2.5; % fused lasso penalty weight

%% set algorithm options (this block doesn't need to be touched)
% termination criterion
options.termin.maxiter = 1000;   % <- maximum number of iterations
options.termin.tol = 1e-6;      % <- relative change in the primal variable
options.termin.progress = 500;   % <- display "progress" (every k iterations...set to inf to disable)
options.termin.silence = false; % <- display termination condition

%=========================================================================%
% stuffs needed for fused lasso
%=========================================================================%
C=tak_adjmat2incmat(adjmat);

%-------------------------------------------------------------------------%
% spectral norm: normest by far the fastest...(but not the most accurate)
%-------------------------------------------------------------------------%
options.sigma=1; % CPPD parameter (sigma*tau L^2 < 1 must be satisfied)
F = [options.lambda*speye(p);options.gamma*C];
options.L=svds(F,1);
options.tau=1/(options.L^2 * options.sigma);
options.tau = options.tau - options.tau/100; % <- safeguard (sig*tau*L^2 < 1...strict equality)
% toc
% return

tic
% options.rho=1;output=tak_admm_enet_regr(Xtr,ytr,options,wtrue);
% output=tak_apgm_flas_regr(Xtr,ytr,options,C,wtrue);
output=tak_cppd_flas_regr(Xtr,ytr,options,C,wtrue);
toc
% fval1=output.fval(1)
% fval_end=output.fval(end)
% [norm(ytr)^2,norm(ytr-Xtr*output.w)^2]
%%
west=output.w;
purge
% nnz(west)
figure,imexptl
subplot(131),tplot(log10(output.fval(2:end))), title('log10(function value)')
subplot(132),tplot(log10(output.wdist)), title('log10(||wtrue-west||)')
subplot(133),tplot(log10(output.rel_changevec)), title('log10(wnew-wold)')

figure,imexptl
subplot(121),imcov(tak_dvecinv(wtrue,0)); title(['wtrue (nnz=',num2str(nnz(wtrue)),')'])
subplot(122),imcov(tak_dvecinv(west,0));  title(['west (nnz=',num2str(nnz(west)),')'])
% subplot(133),imcov(tak_dvecinv(west,0)~=0);  title(['west (support)'])

figure,imexptl
subplot(141),tplot(wtrue)
subplot(142),tplot(west)
subplot(143),tplot2(wtrue,west), legend('wtrue','west'),title('')
subplot(144),tplot(abs(wtrue-west)),title('|wtrue-west|')
% tplottl(abs(w+output.v))
drawnow
%% testing
mse_test = norm(yts_clean - Xts*west)

figure,imexpl
% subplot(411),tplot(yts_clean)
% subplot(412),tplot(Xts*west)
subplot(121),tstem2(yts_clean,Xts*west), legend('yts_clean','yts_est'),title('')
subplot(122),tstem(abs(yts_clean-Xts*west)),title('|yts_clean-yts_est|')
% tplottl(abs(w+output.v))
drawnow