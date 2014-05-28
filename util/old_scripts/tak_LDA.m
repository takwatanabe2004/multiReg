function [ypred,summary] = tak_LDA(Xtrain, ytrain, Xtest, ytest)
% following notation from bickel and levina 2004
% (05/06/2013)
%%
n = length(ytrain);

% n1: +1 group
% n0: -1 group
n1 = length(ytrain==+1);
n0 = length(ytrain==-1);

% Xp = +1 group, Xn = -1 group
Xtrain1 = Xtrain(ytrain==+1,:);
Xtrain0 = Xtrain(ytrain==-1,:);

mu1 = mean(Xtrain1)';
mu0 = mean(Xtrain0)';

DELTA = mu1-mu0;
mu    = (mu0+mu1)/2;

% class probability
pi1 = n1/n;
pi0 = n0/n;

% pooled covariance matrix
SIGMA = (n1*cov(Xtrain1) +  n0*cov(Xtrain0))/(n1+n0);
%% prediction
ntest = length(ytest);
MUrep = repmat(mu,[1,ntest]);

ypred = sign(DELTA'*(SIGMA\(Xtest'-MUrep))) + log(pi1/pi0);
ypred = ypred(:);

% accuracy = sum(ypred==ytest)/ntest;

if nargout > 1
    summary = tak_binary_classification_summary(ypred,ytest);
end

% tplot(DELTA)