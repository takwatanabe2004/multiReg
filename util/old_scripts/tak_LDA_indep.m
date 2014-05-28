function [ypred,summary] = tak_LDA_indep(Xtrain, ytrain, Xtest, ytest)
% following notation from bickel and levina 2004
% independence rule
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

% class probability
pi1 = n1/n;
pi0 = n0/n;

% pooled covariance matrix
SIGMA = (n1*var(Xtrain1) +  n0*var(Xtrain0))/(n1+n0);

%==================================================================================
% remove locations of zero variance
% (this should be fixed once the 'cleansing' stream removes the colinear seeds)
%==================================================================================
idxVar0 = SIGMA==0;
SIGMA(idxVar0) = [];
Xtrain1(:,idxVar0) = [];
Xtrain0(:,idxVar0) = [];
Xtest(:,idxVar0) = [];

p = size(Xtrain1,2);
SIGMAINV = spdiags(SIGMA(:).^-1,0,p,p);

%==================================================================================
% class separation
%==================================================================================
mu1 = mean(Xtrain1)';
mu0 = mean(Xtrain0)';

DELTA = mu1-mu0;
mu    = (mu0+mu1)/2;
%% prediction
ntest = length(ytest);
MUrep = repmat(mu,[1,ntest]);

ypred = SIGN(DELTA'*(SIGMAINV*(Xtest'-MUrep))) + log(pi1/pi0);
ypred = ypred(:);

% accuracy = sum(ypred==ytest)/ntest;
if nargout > 1
    summary = tak_binary_classification_summary(ypred,ytest);
end

% tplot(DELTA)