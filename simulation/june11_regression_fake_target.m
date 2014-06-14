%% june11_regression_fake_target
% (06/11/2014)
%=========================================================================%
% apply regression on fake data made from june11_make_fake_target.m
%=========================================================================%
clear
purge

grid='WashU'; % {'Grid326','Grid1068','WashU'}

load([get_rootdir,'/data_local/designMatrix_FC_',grid,'.mat'],'X')
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',...
    grid,'_dilated5mm.mat'], 'roiLabel', 'yeoLabels')
[n,p]=size(X);
%% load fake weight vector
snr = 0;
load([fileparts(mfilename('fullpath')),...
    '/fake_target_',grid,'_snr',num2str(snr),'.mat'],'y','wtrue')
% tplott(y)
% tplott(wtrue)
% return

%-------------------------------------------------------------------------%
% (optionally create new design matrix)
%-------------------------------------------------------------------------%
n=333;
X=randn(n,p);
X=zscore(X);
snr=0;
mask_supp = wtrue~=0;
wtrue(mask_supp)=wtrue(mask_supp)+5;
y=X*wtrue + snr*randn(n,1);
%% some plot options
cbarOption={'fontsize',22','fontweight','b','ytick',[-.66,0,.66],...
    'YTickLabel',{' <0',' =0',' >0'},'TickLength',[0 0]};
textOption1={'fontweight','b','fontsize',9};
lineOption = {'color','k','linewidth',0.5};
lwidth_deg=2.5;

% get yeo info
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',grid,'_dilated5mm.mat'])

% circularly shift 1 indices (so "unlabeled" is at the final label index)
roiLabel=roiLabel-1;
roiLabel(roiLabel==-1)=12;
yeoLabels=circshift(yeoLabels,-1);

% get sorting indices for visualization
[idxsrt,labelCount] = tak_get_yeo_sort(roiLabel);
%%
%=========================================================================%
% set options
%=========================================================================%
% options.lambda=2^12;  % L1 penalty weight
% options.gamma =2^8; % fused lasso penalty weight

options.lambda=2^0.25;  % L1 penalty weight
options.gamma =2^5; % fused lasso penalty weight
%% set algorithm options (this block doesn't need to be touched)
% termination criterion
options.termin.maxiter = 500;   % <- maximum number of iterations
options.termin.tol = 5e-4;      % <- relative change in the primal variable
options.termin.progress = 20;   % <- display "progress" (every k iterations...set to inf to disable)
options.termin.silence = false; % <- display termination condition
options.fval = true; % <- keep track of function values (may slow down algorithm)

%=========================================================================%
% stuffs needed for fused lasso & CPPD alg
%=========================================================================%
% %-------------------------------------------------------------------------%
% % spectral norm for step size
% %-------------------------------------------------------------------------%
% load([get_rootdir,'/data_local/graphinfo/graph_info_',grid,'.mat'],'C')
% % C=tak_adjmat2incmat(adjmat);
% 
% options.sigma=1; % CPPD parameter (sigma*tau L^2 < 1 must be satisfied)
% F = [options.lambda*speye(p);options.gamma*C];
% % tic
% % options.L=svds(F,1)
% options.L=sqrt(eigs(F'*F,1));
% % toc
% options.tau=1/(options.L^2 * options.sigma)
% options.tau = options.tau - options.tau/100; % <- safeguard (sig*tau*L^2 < 1...strict equality)
% % toc
% % return


tic
options.rho=1;output=tak_admm_enet_regr(X,y,options,wtrue);
% output=tak_apgm_flas_regr(Xtr,ytr,options,C,wtrue);
% output=tak_cppd_flas_regr(X,y,options,C,wtrue);
toc
% fval1=output.fval(1)
% fval_end=output.fval(end)
% [norm(ytr)^2,norm(ytr-Xtr*output.w)^2]
%%
west=output.w;
west_mat = tak_dvecinv(west,0);
west_mat_srt = west_mat(idxsrt,idxsrt);

figure,imexp
subplot(121)
    imcov(west_mat_srt),axis off,
    tmp=max(abs(caxis));    caxis([-tmp,tmp]/.5)
    colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
subplot(122)
    wtrue_mat=tak_dvecinv(wtrue,0);
    imcov(wtrue_mat(idxsrt,idxsrt)),axis off,
    tmp=max(abs(caxis));    caxis([-tmp,tmp]/1)
    colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
figure,imexpl
    Xty_mat=tak_dvecinv(X'*y,0);
    imcov(Xty_mat(idxsrt,idxsrt)),axis off,
    tmp=max(abs(caxis));    caxis([-tmp,tmp]/1)
    colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)


figure,imexp
% subplot(311),tstem(y),title('y')
% subplot(312),tstem(X*west),title('yest')
subplot(211),tstem2(y,X*west),legend('y','yest'),YLIM=ylim;
subplot(212),tstem(abs(y-X*west)),title('|y-yest|'),ylim(YLIM)
drawnow


% mse_test = norm(y - X*west)
corr_test = corr(y,X*west)

% nnz(west)
figure,imexp
subplot(241),tplot(log10(output.fval(2:end))), title('log10(function value)')
subplot(242),tplot(log10(output.wdist)), title('log10(||wtrue-west||)')
subplot(243),tplot(log10(output.rel_changevec)), title('log10(wnew-wold)')

subplot(245),tplot(wtrue)
subplot(246),tplot(west)
subplot(247),tplot2(wtrue,west), legend('wtrue','west'),title('')
subplot(248),tplot(abs(wtrue-west)),title('|wtrue-west|')
% tplottl(abs(w+output.v))
% drawnow

%%
% loss_val = norm(age-X*west)^2/2
% lasso_val = options.lambda * norm(west,1)
% flass_val = options.gamma * norm(C*west,1)
% enet_val  = options.gamma * norm(west)^2