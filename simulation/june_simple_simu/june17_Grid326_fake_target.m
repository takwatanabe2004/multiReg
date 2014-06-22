%% june17_Grid326_fake_target
% (06/17/2014)
%=========================================================================%
% - 2-d fused lasso trial
%=========================================================================%
%%
clear all;
purge
randn('state',0)

GRID='Grid326'; % {'Grid326','Grid1068','WashU'}

load([get_rootdir,'/data_local/graphinfo/graph_info_',GRID,'.mat'],'C')

load([get_rootdir,'/data_local/designMatrix_FC_',GRID,'.mat'],'X')
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',...
    GRID,'_dilated5mm.mat'], 'roiLabel', 'yeoLabels')
[~,p]=size(X);
%% load fake weight vector
sig = 0;
load([get_rootdir,'/simulation',...
    '/fake_target_',GRID,'_snr',num2str(sig),'.mat'],'y','wtrue')
% tplott(y)
% tplott(wtrue)
% return

%-------------------------------------------------------------------------%
% (optionally create new design matrix)
%-------------------------------------------------------------------------%
n=2000;
X=randn(n,p);
X=zscore(X);
sig=0;
mask_supp = wtrue~=0;
wtrue(mask_supp)=wtrue(mask_supp)+5;
y=X*wtrue + sig*randn(n,1);

ntest=500;
Xtest=randn(ntest,p);
ytest=Xtest*wtrue;
% return
%% set options for admm
options_EN.rho=1; % <- AL parameter
options_EN.maxiter = 500;      % <- maximum number of iterations
options_EN.tol = 1e-3;         % <- relative change in the primal variable
options_EN.progress = inf;     % <- display "progress" (every k iterations...set to inf to disable)
options_EN.silence = false;    % <- display termination condition
options_EN.funcval = false;    % <- track function values (may slow alg.)

options_FL = options_EN;
%% some plot options
cbarOption={'fontsize',22','fontweight','b','ytick',[-.66,0,.66],...
    'YTickLabel',{' <0',' =0',' >0'},'TickLength',[0 0]};
textOption1={'fontweight','b','fontsize',9};
lineOption = {'color','k','linewidth',0.5};
lwidth_deg=2.5;

% get yeo info
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',GRID,'_dilated5mm.mat'])

% circularly shift 1 indices (so "unlabeled" is at the final label index)
roiLabel=roiLabel-1;
roiLabel(roiLabel==-1)=12;
yeoLabels=circshift(yeoLabels,-1);

% get sorting indices for visualization
[idxsrt,labelCount] = tak_get_yeo_sort(roiLabel);
%% estimate weight vector
%-------------------------------------------------------------------------%
% ridge regression
%-------------------------------------------------------------------------%
% w_RR = tak_ridge_regression(X,y,.05);
%%
%-------------------------------------------------------------------------%
% enet
%-------------------------------------------------------------------------%
lam1 = 1;  % L1 penalty weight
gam1 = 0; % fused lasso penalty weight
if ~isfield(options_EN, 'K'), 
    tic
    options_EN.K = tak_admm_inv_lemma(X,1/(options_EN.rho+gam1));
    toc
end
% [w_EN, output_EN]=tak_admm_EN_regr(X,y,lam1,gam1,options_EN,wtrue(:));

%
%-------------------------------------------------------------------------%
% fused lasso
%-------------------------------------------------------------------------%
lam2 = .1;  % L1 penalty weight
gam2 = 12; % fused lasso penalty weight
if ~isfield(options_FL, 'K'), 
    tic
    options_FL.K = tak_admm_inv_lemma(X,1/options_FL.rho);
    toc
end
[w_FL,output_FL]=tak_admm_FL_regr_pcg(X,y,lam2,gam2,options_FL,C,[],wtrue(:));

% figure,imexp
% subplot(241),tplot(log10(output_FL.fval(2:end))), title('log10(function value)')
% subplot(242),tplot(log10(output_FL.wdist)), title('log10(||wtrue-west||)')
% subplot(243),tplot(log10(output_FL.rel_changevec)), title('log10(wnew-wold)')
% 
% subplot(245),tplot(wtrue)
% subplot(246),tplot(w_FL)
% subplot(247),tplot2(wtrue,w_FL), legend('wtrue','west'),title('')
% subplot(248),tplot(abs(wtrue-w_FL)),title('|wtrue-west|')
%%% compare methods
% purge
wtrue_mat = tak_dvecinv(wtrue, 0); wtrue_mat=wtrue_mat(idxsrt,idxsrt);

cscale=.5;
figure,imexpb
subplot(141)
    imcov(wtrue_mat),axis off,
    tmp=max(abs(caxis));    caxis([-tmp,tmp]*cscale)
    colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
if exist('w_RR','var')
    subplot(142)
    w_RR_mat = tak_dvecinv(w_RR, 0); w_RR_mat=w_RR_mat(idxsrt,idxsrt);
    imcov(w_RR_mat),axis off,
    tmp=max(abs(caxis));    caxis([-tmp,tmp]*cscale)
    colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
end
if exist('w_EN','var')
    subplot(143)
    w_EN_mat = tak_dvecinv(w_EN, 0); w_EN_mat=w_EN_mat(idxsrt,idxsrt);
    imcov(w_EN_mat),axis off,
    tmp=max(abs(caxis));    caxis([-tmp,tmp]*cscale)
    colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
end
if exist('w_FL','var')
    subplot(144)
    w_FL_mat = tak_dvecinv(w_FL, 0); w_FL_mat=w_FL_mat(idxsrt,idxsrt);
    imcov(w_FL_mat),axis off,
    tmp=max(abs(caxis));    caxis([-tmp,tmp]*cscale)
    colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
end
lwid=2;
% figure,imexpb
% % subplot(221),imcov(wtrue),caxis(CAXIS); colorbar('location','northoutside')
% stem(ytest,'linewidth',lwid),xlim([1,length(ytest)]), hold on
% stem(Xtest*w_RR,'linewidth',lwid,'color','k')
% stem(Xtest*w_EN,'linewidth',lwid,'color','g')
% stem(Xtest*w_FL,'linewidth',lwid,'color','r')
% legend('true','RR','EN','FL'), ylim([min(ytest),max(ytest)]),grid on
% drawnow

if exist('w_RR','var'), mse_RR = norm(ytest - Xtest*w_RR),end;
if exist('w_EN','var'), mse_EN = norm(ytest - Xtest*w_EN), end;
if exist('w_FL','var'), mse_FL = norm(ytest - Xtest*w_FL), end;
if exist('w_RR','var'), corr_RR = corr(ytest, Xtest*w_RR), end;
if exist('w_EN','var'), corr_EN = corr(ytest, Xtest*w_EN), end;
if exist('w_FL','var'), corr_FL = corr(ytest, Xtest*w_FL), end;