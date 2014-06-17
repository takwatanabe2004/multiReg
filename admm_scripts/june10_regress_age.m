%% june10_regress_age
% (06/10/2014)
%=========================================================================%
% study feature matrix
%=========================================================================%
clear
purge

grid='Grid326'; % {'Grid326','Grid1068','WashU'}

load([get_rootdir,'/data_local/designMatrix_FC_',grid,'.mat'],'X','age')
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',grid,'_dilated5mm.mat'])

[n,p]=size(X);
% tplot(sort(age))
% figure,hist(age,[min(age):max(age)])
% return

%=========================================================================%
% stuffs needed for fused lasso
%=========================================================================%
load([get_rootdir,'/data_local/graphinfo/graph_info_',grid,'.mat'],'C')
% C=tak_adjmat2incmat(adjmat);
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
[idxsrt,labelCount] = tak_get_yeo_sort(roiLabel);
%%
%=========================================================================%
% set options
%=========================================================================%
% options.lambda=0;  % L1 penalty weight
% options.gamma =5; % fused lasso penalty weight
options.lambda=2^4.5;  % L1 penalty weight
options.gamma =2^-4.5; % fused lasso penalty weight
%% set algorithm options (this block doesn't need to be touched)
% termination criterion
options.termin.maxiter = 500;   % <- maximum number of iterations
options.termin.tol = 5e-4;      % <- relative change in the primal variable
options.termin.progress = 20;   % <- display "progress" (every k iterations...set to inf to disable)
options.termin.silence = false; % <- display termination condition
options.fval = true; % <- keep track of function values (may slow down algorithm)

%-------------------------------------------------------------------------%
% spectral norm for step size
%-------------------------------------------------------------------------%
options.sigma=1; % CPPD parameter (sigma*tau L^2 < 1 must be satisfied)
F = [options.lambda*speye(p);options.gamma*C];
% tic
% options.L=svds(F,1)
options.L=sqrt(eigs(F'*F,1));
% toc
options.tau=1/(options.L^2 * options.sigma)
options.tau = options.tau - options.tau/100; % <- safeguard (sig*tau*L^2 < 1...strict equality)
toc
return


tic
options.rho=1;output=tak_admm_enet_regr(X,age,options);
% output=tak_apgm_flas_regr(Xtr,ytr,options,C,wtrue);
% output=tak_cppd_flas_regr(X,age,options,C);
toc
% fval1=output.fval(1)
% fval_end=output.fval(end)
% [norm(ytr)^2,norm(ytr-Xtr*output.w)^2]
%%
west=output.w;
west_mat = tak_dvecinv(west,0);
west_mat_srt = west_mat(idxsrt,idxsrt);

purge
% nnz(west)
figure,imexptl
subplot(131),tplot(log10(output.fval(2:end))), title('log10(function value)')
subplot(132),tplot(log10(output.rel_changevec)), title('log10(wnew-wold)')
subplot(133),tplot(west)

% figure,imexptl
% subplot(121),imcov(west_mat);  title(['west (nnz=',num2str(nnz(west)),')'])
%     tmp=max(abs(caxis));    caxis([-tmp,tmp]/15)
% subplot(122),imcov(west_mat_srt);  title(['west (nnz=',num2str(nnz(west)),')'])
%     tmp=max(abs(caxis));    caxis([-tmp,tmp]/15)
figure,imexpl
    imcov(west_mat_srt),axis off,
    tmp=max(abs(caxis));    caxis([-tmp,tmp]/5)
    colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)

% mse_test = norm(age - X*west)
corr_test = corr(age,X*west)

figure,imexpt
% subplot(311),tstem(age), title('age')
% subplot(312),tstem(X*west), title('age_est')
subplot(211),tstem2(age,X*west), legend('age','ageEst'),title('')
subplot(212),tstem(abs(age-X*west)),title('|age-age_est|')
% tplottl(abs(w+output.v))
drawnow
%%
% loss_val = norm(age-X*west)^2/2
% lasso_val = options.lambda * norm(west,1)
% flass_val = options.gamma * norm(C*west,1)
% enet_val  = options.gamma * norm(west)^2