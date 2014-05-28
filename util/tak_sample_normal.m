function X = tak_sample_normal(nsamp, covmat)
% X = tak_sample_normal(nsamp, covmat)
%|--------------------------------------------------------------------------------------|%
%| Description:
%|      Sample from the normal distribution parameterized by the 
%|      covariance matrix.
%|--------------------------------------------------------------------------------------|%
%| Input: 
%|      nsamp = number of samples to draw
%|       covmat = covariance matrix
%|--------------------------------------------------------------------------------------|%
%| Output: 
%|          X = (nsamp x p) data matrix
%|--------------------------------------------------------------------------------------|%
%| Created 8/01/2012
%| 04/14/2013 -> X is now (nsamp x p) matrix
%|--------------------------------------------------------------------------------------|%
%%
% warning('4/14/2013 -> X is now an (n x p) matrix')
p = size(covmat,1);

%| compute cholesky
L = chol(covmat,'lower');

z = randn(p,nsamp);

%| draw realizations
X = (L*z)';