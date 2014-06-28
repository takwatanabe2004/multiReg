%% study_PCA_rcorr.m
% (05/12/2014)
%=========================================================================%
% - Visualize few of the PC's in the FC in the connectome space
%=========================================================================%
%%
clear
purge

%% load data
load rcorr_design_censor.mat X y

[n,p]=size(X);

% scale to 0 mean, unit variance?
flag_scale = false;
if flag_scale
    X=zscore(X);
end
%% do PCA
%=========================================================================%
% Do PCA: 
%   Since rank(X) = n, project X into the n-d variance maximizing subspace.
%   To get the first k principal coordinates, do X(:,1:k), where k<=n.
%
%   Variance normalization is not done, since all features are in the same scale
%-------------------------------------------------------------------------%
% X  (in): (n x p) data in original feature space
% X (out): (n x n) data projected into the n-d variance maximizing subspace
% SVs: (p x 1) vector of singular values.  
%       singular_values.^2/n is the principal component values 
%       (ie, eigenvalues of the empirical covariance)
%=========================================================================%
tic
X = tak_center_features(X);
[U,SVs,V] = svd(X, 'econ');
tPCA = toc
SVs = diag(SVs);

% encode features into the variance maximizing subspace
X = X*V;
%% visualize PC loadings
purge
%==========================================================================
% yeo network coordinate info (for visualization)
%==========================================================================
load yeo_info347_dilated5mm.mat yeoLabels roiMNI roiLabel
[idx_sort,labelCount] = tak_get_yeo_sort(roiLabel);

% circularly shift 1 indices
roiLabel=roiLabel-1;
roiLabel(roiLabel==-1)=12;

wmat1=tak_dvecinv(V(:,1),0);
wmatsort1=wmat1(idx_sort,idx_sort);

wmat2=tak_dvecinv(V(:,2),0);
wmatsort2=wmat2(idx_sort,idx_sort);

wmat3=tak_dvecinv(V(:,3),0);
wmatsort3=wmat3(idx_sort,idx_sort);

wmat4=tak_dvecinv(V(:,4),0);
wmatsort4=wmat4(idx_sort,idx_sort);

wmat120=tak_dvecinv(V(:,120),0);
wmatsort120=wmat120(idx_sort,idx_sort);

wmat119=tak_dvecinv(V(:,119),0);
wmatsort119=wmat119(idx_sort,idx_sort);

%-------------------------------------------------------------------------%
% some figure options
%-------------------------------------------------------------------------%
% text/line option for the line-boxes for the network partitioning
textOption1={'fontweight','b','fontsize',12};
lineOption = {'color','k','linewidth',0.5};

%=========================================================================%
% show figures
%=========================================================================%
figure,imexp
subplot(121),imcov(wmatsort1),axis off
tak_local_linegroups3(gcf,labelCount,textOption1,yeoLabels,lineOption),hold on
subplot(122),imcov(wmatsort2),axis off
tak_local_linegroups3(gcf,labelCount,textOption1,yeoLabels,lineOption),hold on

figure,imexp
subplot(121),imcov(wmatsort3),axis off
tak_local_linegroups3(gcf,labelCount,textOption1,yeoLabels,lineOption),hold on
subplot(122),imcov(wmatsort4),axis off
tak_local_linegroups3(gcf,labelCount,textOption1,yeoLabels,lineOption),hold on

figure,imexp
subplot(121),imcov(wmatsort119),axis off
tak_local_linegroups3(gcf,labelCount,textOption1,yeoLabels,lineOption),hold on
subplot(122),imcov(wmatsort120),axis off
tak_local_linegroups3(gcf,labelCount,textOption1,yeoLabels,lineOption),hold on
