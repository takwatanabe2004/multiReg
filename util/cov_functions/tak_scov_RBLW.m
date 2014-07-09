function [scovRBLW, shrinkage] = tak_scov_RBLW(X)
%| X = (n x p) design matrix
%| Ref: 2010 Chen, Hero etal "Shrinkage algorithms for mmse covariance estimation"
%| 12/17/2012
%%
%| sample covariance normalized by n-1
scov = cov(X,1);

[n,p] = size(X);

%| compute numerator term of the shrinkage coefficient
num = (n-2)/n * trace(scov^2) + trace(scov)^2;

%| denominator term
den = (n+2) * (trace(scov^2) - trace(scov)^2/p);

%| the shrinkage coefficient
shrinkage = min(num/den,1);

%| the shrinkage target 
target = trace(scov)/p*eye(p);
scovRBLW = (1-shrinkage)*scov + shrinkage*target;