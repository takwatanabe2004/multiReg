function [scovLW, shrinkage] = tak_scov_LW(X)
%| X = (n x p) design matrix
%| Ref: 2010 Chen, Hero etal "Shrinkage algorithms for mmse covariance estimation"
%| 12/17/2012
%%
%| sample covariance normalized by n
scov = cov(X,1);

[n,p] = size(X);

%| compute numerator term of the shrinkage coefficient
num = 0;
for i = 1:n
    xi = X(i,:);
    num = num + norm(xi*xi' - scov,'fro')^2;
end

%| denominator term
den = n^2 * (trace(scov^2) - trace(scov)^2/p);

%| the shrinkage coefficient
shrinkage = min(num/den,1);

%| the shrinkage target 
target = trace(scov)/p*eye(p);
scovLW = (1-shrinkage)*scov + shrinkage*target;