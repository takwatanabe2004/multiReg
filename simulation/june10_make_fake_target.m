%% june10_make_fake_target.m
% (06/10/2014)
%=========================================================================%
% Make artificial target data (y) for regression
%=========================================================================%
clear
purge

grid='Grid326'; % {'Grid326','Grid1068','WashU'}

load([get_rootdir,'/data_local/designMatrix_FC_',grid,'.mat'],'X')
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',grid,'_dilated5mm.mat'])

% [age,idx_forw]=sort(age);
% [~,idx_back]=sort(idx_forw);
[n,p]=size(X);
%% make fake target
% get yeo info
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',grid,'_dilated5mm.mat'],...
    'roiLabel','yeoLabels')
d=length(roiLabel);
%% some plot options
textOption1={'fontweight','b','fontsize',12};
lineOption = {'color','k','linewidth',0.5};
lwidth_deg=2.5;

% circularly shift 1 indices (so "unlabeled" is at the final label index)
roiLabel=roiLabel-1;
roiLabel(roiLabel==-1)=12;
yeoLabels=circshift(yeoLabels,-1);
[idxsrt,labelCount] = tak_get_yeo_sort(roiLabel);
[~,idxsrt_back] = sort(idxsrt);
% idxsrt(idxsrt_back)
% return
%%
roiLabel_srt = roiLabel(idxsrt);
mask_mat = false(d);

% mask7=(roiLabel_srt==7);
% mask6=(roiLabel_srt==6);
% % mask = (roiLabel_srt==6) | (roiLabel_srt==7);

for i = 1:d
    f1i = roiLabel_srt(i)==1-1;
    f2i = roiLabel_srt(i)==2-1;
    f6i = roiLabel_srt(i)==6-1;
    f7i = roiLabel_srt(i)==7-1;
    for j=1+i:d
        f1j = roiLabel_srt(j)==1-1;
        f2j = roiLabel_srt(j)==2-1;
        f3j = roiLabel_srt(j)==3-1;
        f6j = roiLabel_srt(j)==6-1;
        f7j = roiLabel_srt(j)==7-1;
        f12j = roiLabel_srt(j)==12-1;
        if (f1i & f2j) | (f2i & f3j)  | (f2i & f12j) | (f7i&f7j)
            mask_mat(i,j) = true;
        end
    end
end
mask_mat=mask_mat';        
mask_vec=tak_dvec(mask_mat);
% imedgel(mask_mat+eye(d),0)
%     tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
wtrue_srt = zeros(p,1);
wtrue_srt(mask_vec)=randn(sum(mask_vec),1) + 10000;
wtrue_srt_mat = tak_dvecinv(wtrue_srt,0);
% imcovvl(wtrue_srt_mat),  axis off,    colorbar('location','northoutside')
%     tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
wtrue_mat = wtrue_srt_mat(idxsrt_back,idxsrt_back);
wtrue = tak_dvec(wtrue_mat);
% imcovvl(wtrue_mat)
% imcovvl(wtrue_mat(idxsrt,idxsrt))
% return
%% make fake data
snr=0;
X = randn(n,p);
y = X*wtrue + snr*randn(n,1);

% return
% X=X(idx,:);
% X0=X;
% meanX=mean(X0,1);
% X=zscore(X);

% lambda=55555e3
lambda=.1e10
%%
% compute matrix K...an expression frequently encountered during admm
% when using the inversion lemma
%   K = X'*inv( I + c*X*X') 
% X: (n x p) design matrix
Xty = X'*y;
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

mse_test = norm(y - X*west)/norm(y)
corr_test = corr(y,X*west)

figure,imexpr
% subplot(311),tstem(age), title('age')
% subplot(312),tstem(X*west), title('age_est')
subplot(211),tstem2(y,X*west), legend('age','ageEst'),title('')
subplot(212),tstem(abs(y-X*west)),title('|age-age_est|')
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

imcovvl(wtrue_mat(idxsrt,idxsrt))
    colorbar('location','northoutside')
    tmp=max(abs(caxis));    caxis([-tmp,tmp]/1)
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
%%
% loss_val = norm(age-X*west)^2/2
% lasso_val = options.lambda * norm(west,1)
% flass_val = options.gamma * norm(C*west,1)
% enet_val  = options.gamma * norm(west)^2