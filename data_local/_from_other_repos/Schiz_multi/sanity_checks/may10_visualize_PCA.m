%% may10_visualize_PCA
% (05/10/2014)
%=========================================================================%
% - visualize pca loadings
%=========================================================================%
clear
purge
%% load data
load sMRI_design_censor.mat X y

[n,p]=size(X);

% scale to 0 mean, unit variance?
flag_scale = true;
if flag_scale
    X=zscore(X);
end
%% do pca
tic
Xctr = tak_center_features(X);
[U,SVs,V] = svd(Xctr, 'econ');
tPCA = toc
SVs = diag(SVs);

% encode features into the variance maximizing subspace
X = Xctr*V;
%%
% visualize PC loadings
purge
load([get_rootdir, '/data_local/A_matrix_fullFOV.mat'], 'A_full','nx','ny','nz')
nii =tak_downsample_nii(load_nii([get_rootdir,'/data_local/m0wrp1mprage.nii']));

% first component
v1 = reshape(A_full*V(:,1), [nx,ny,nz]);
v2 = reshape(A_full*V(:,2), [nx,ny,nz]);
v3 = reshape(A_full*V(:,3), [nx,ny,nz]);
% tak_gui_show_slices(v1)
% tak_gui_show_slices(v2)
% tak_gui_show_slices(v3)

nii_v1=nii;
nii_v2=nii;
nii_v3=nii;
nii_v1.img=v1;
nii_v2.img=v2;
nii_v3.img=v3;
view_nii(nii_v1);
view_nii(nii_v2);
view_nii(nii_v3);

vv = reshape(A_full*V(:,111), [nx,ny,nz]);
nii_vv=nii;
nii_vv.img=vv;
view_nii(nii_vv);