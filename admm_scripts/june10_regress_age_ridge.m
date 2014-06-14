%% june10_regress_age_ridge
% (06/10/2014)
%=========================================================================%
% Regress age using ridge-regression
%=========================================================================%
clear
purge

grid='Grid326'; % {'Grid326','Grid1068','WashU'}

load([get_rootdir,'/data_local/designMatrix_FC_',grid,'.mat'],'X','age')
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',grid,'_dilated5mm.mat'])

[age,idx]=sort(age);
[n,p]=size(X);

% X=X(idx,:);
% X0=X;
% meanX=mean(X0,1);
% X=zscore(X);

% lambda=55555e3
lambda=1e-4
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
% compute matrix K...an expression frequently encountered during admm
% when using the inversion lemma
%   K = X'*inv( I + c*X*X') 
% X: (n x p) design matrix
Xty = X'*age;
TMP= eye(n) + (X*X')/lambda;
K = (TMP\X)';
west = Xty/lambda - 1/lambda^2*K*(X*Xty);
%%
west_mat = tak_dvecinv(west,0);
west_mat_srt = west_mat(idxsrt,idxsrt);

purge
% nnz(west)
% figure,imexptl
% subplot(131),tplot(log10(output.fval(2:end))), title('log10(function value)')
% subplot(132),tplot(log10(output.rel_changevec)), title('log10(wnew-wold)')
% subplot(133),tplot(west)

% figure,imexptl
% subplot(121),imcov(west_mat);  title(['west (nnz=',num2str(nnz(west)),')'])
%     tmp=max(abs(caxis));    caxis([-tmp,tmp]/15)
% subplot(122),imcov(west_mat_srt);  title(['west (nnz=',num2str(nnz(west)),')'])
%     tmp=max(abs(caxis));    caxis([-tmp,tmp]/15)
figure,imexpr
    imcov(west_mat_srt),axis off,
    tmp=max(abs(caxis));    caxis([-tmp,tmp]/1)
    colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)

mse_test = norm(age - X*west)/norm(age)
corr_test = corr(age,X*west)

figure,imexpr
% subplot(311),tstem(age), title('age')
% subplot(312),tstem(X*west), title('age_est')
subplot(211),tstem2(age,X*west), legend('age','ageEst'),title('')
subplot(212),tstem(abs(age-X*west)),title('|age-age_est|')
% tplottl(abs(w+output.v))
drawnow
%%
Xty_mat = tak_dvecinv(Xty,0);
Xty_mat_srt=Xty_mat(idxsrt,idxsrt);
% imcovvl(Xty_mat)
imcovvl(Xty_mat_srt),axis off
    colorbar('location','northoutside')
    tmp=max(abs(caxis));    caxis([-tmp,tmp]/1)
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
%%
% loss_val = norm(age-X*west)^2/2
% lasso_val = options.lambda * norm(west,1)
% flass_val = options.gamma * norm(C*west,1)
% enet_val  = options.gamma * norm(west)^2