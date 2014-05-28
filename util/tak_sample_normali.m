function X = tak_sample_normali(nsamp, icov)
% X = tak_sample_normali(nsamp, icov)
%|--------------------------------------------------------------------------------------|%
%| Description:
%|      Sample from the normal distribution parameterized by the 
%|      inverse covariance matrix.
%|--------------------------------------------------------------------------------------|%
%| Input: 
%|      nsamp = number of samples to draw
%|       icov = inverse covariance matrix
%|--------------------------------------------------------------------------------------|%
%| Output: 
%|          X = (nsamp x p) data matrix
%|--------------------------------------------------------------------------------------|%
%| Created 8/01/2012
%| 04/14/2013 -> X is now (nsamp x p) matrix
%|--------------------------------------------------------------------------------------|%
% warning('4/14/2013 -> X is now an (n x p) matrix')
p = size(icov,1);

%| compute cholesky
% 06/10/2013 -> just noticed I used letter L to denote the upper cholesky factor...
%               which was probably a dumb idea...implementation is correct, but 
%               confused myself today....
L = chol(icov,'upper');

z = randn(p,nsamp);

%| draw realizations
X = (L\z)';