function summary = tak_binary_classification_summary_extended(ypred, ytrue)
%| summary = tak_binary_classification_summary_extended(ypred, ytrue)
%|--------------------------------------------------------------------------------|
%| 'description'
%|--------------------------------------------------------------------------------|
%| INPUT
%|   in1: 'description'
%|   in2: 'description'
%|--------------------------------------------------------------------------------|
%| OUTPUT
%|  summary:
%|--------------------------------------------------------------------------------|
%| Example:
%
%|--------------------------------------------------------------------------------|
%| (12/25/2013)
%| - returns more info than the original one
%%
if length(ytrue)~= length(ypred)
    error('dimension mismatch')
end

% (09/19/2013)
ytrue=ytrue(:);
ypred=ypred(:);

n = size(ytrue,1);

%==================================================================================
% check if the labels are all +1 or -1
%==================================================================================
if n~=( sum(ytrue==1)+sum(ytrue==-1) ) || n~=( sum(ypred==1)+sum(ypred==-1) )
    error(['The entries of ''ytrue'' and ''ypred'' must all be +1 or -1...',...
             'return only the overall accuracy'])
end

% overall accuracy
summary.accuracy  = sum(ytrue==ypred)/n;

% 'p' for positive, 'n' for negative
idxp = (ytrue == +1);
idxn = (ytrue == -1);

%==================================================================================
% TP = true positives
% TN = true negatives
% FP = false positives
% FN = false negatives
%==================================================================================
TP = sum(ytrue(idxp)==ypred(idxp));
TN = sum(ytrue(idxn)==ypred(idxn));
FP = sum(ytrue(idxn)~=ypred(idxn));
FN = sum(ytrue(idxp)~=ypred(idxp));

%==================================================================================
% TPR = true positive rate  (aka sensitivity, recall, hit rate, power, detection rate)
% TNR = true negative rate  (aka specificity)
% FPR = false positive rate (aka size, type I error rate)
% FNR = false negative rate (aka miss rate, type II error rate, 1-TPR)
% PPV = positive predictive value (aka precision)
% NPV = negative predictive value
%==================================================================================
summary.TPR = TP/(TP+FN);
summary.TNR = TN/(TN+FP);

summary.FPR = FP/(FP+TN);
summary.FNR = FN/(FN+TP);
summary.PPV = TP/(TP+FP);
summary.NPV = TN/(TN+FN);

summary.TP = TP;
summary.TN = TN;
summary.FP = FP;
summary.FN = FN;