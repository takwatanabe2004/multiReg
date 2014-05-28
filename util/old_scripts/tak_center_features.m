function Xctr = tak_center_features(X)
% Xctr = tak_standardize_features(X)
%|------------------------------------------------------------------------------|%
% Center/demean the (n x p) design matrix X so that all p predictors are mean zero.
%|------------------------------------------------------------------------------|%
% INPUT
%      X: (n x p) design matrix (n realizations, p predictors)
%|------------------------------------------------------------------------------|%
% OUTPUT
%  Xstd: (n x p) centered design matrix
%|------------------------------------------------------------------------------|%
%| 03/11/2013
%| 04/09/2013 -> updated to more efficient version
%%
%| from matlab function 'zscore'
Xctr = bsxfun(@minus,X,mean(X));