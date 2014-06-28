%% save_PCA_MM.m
% (05/12/2014)
%=========================================================================%
% - Save the first n principal loading vectors for the (n x p) 
%   design matrix X (sMRI + FC dataset)
%=========================================================================%
%%
clear
purge

fsave=true;
%% load data

load MM_design_censor.mat X y

% scale to 0 mean, unit variance?
flag_scale = 0;
if flag_scale
    X=zscore(X);
end

%% output save path
if flag_scale
    outPath=[get_rootdir,'/data_local/multimodal/PCA_zscaled.mat'];
else
    outPath=[get_rootdir,'/data_local/multimodal/PCA.mat'];
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
%% save
mFileName=mfilename;
timeStamp=tak_timestamp;

if fsave
    save(outPath,outVars{:})
end