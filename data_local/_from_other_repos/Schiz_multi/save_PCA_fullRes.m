%% save_PCA_fullRes.m
% (05/10/2014)
%=========================================================================%
% - Save the first n principal loading vectors for the (n x p) 
%   design matrix X
%=========================================================================%
%%
clear
purge

fsave=true;
%% load data
load sMRI_design_censor_fullRes.mat X y

[n,p]=size(X);

% scale to 0 mean, unit variance?
flag_scale = false;
if flag_scale
    X=zscore(X);
end

%% output save path
if flag_scale
    outPath=[get_rootdir,'/data_local/PCA_zscaled_fullRes.mat'];
else
    outPath=[get_rootdir,'/data_local/PCA_fullRes.mat'];
end
outVars={'X','SVs','timeStamp'};
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
% purge
% load([get_rootdir, '/data_local/A_matrix_fullFOV.mat'], 'A_full','nx','ny','nz')
% nii =tak_downsample_nii(load_nii([get_rootdir,'/data_local/m0wrp1mprage.nii']));
% 
% % first component
% v1 = reshape(A_full*V(:,1), [nx,ny,nz]);
% v2 = reshape(A_full*V(:,2), [nx,ny,nz]);
% v3 = reshape(A_full*V(:,3), [nx,ny,nz]);
% % tak_gui_show_slices(v1)
% % tak_gui_show_slices(v2)
% % tak_gui_show_slices(v3)
% 
% nii_v1=nii;
% nii_v2=nii;
% nii_v3=nii;
% nii_v1.img=v1;
% nii_v2.img=v2;
% nii_v3.img=v3;
% view_nii(nii_v1);
% view_nii(nii_v2);
% view_nii(nii_v3);
% 
% vv = reshape(A_full*V(:,111), [nx,ny,nz]);
% nii_vv=nii;
% nii_vv.img=vv;
% view_nii(nii_vv)
%% save
mFileName=mfilename;
timeStamp=tak_timestamp;

if fsave
    save(outPath,outVars{:})
end