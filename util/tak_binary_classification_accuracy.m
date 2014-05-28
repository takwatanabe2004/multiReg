function accuracy = tak_binary_classification_accuracy(ypred, ytrue)
% summary = tak_binary_classification_accuracy(ypred, ytrue,option)
% - report accuracy
%-------------------------------------------------------------------------%
% (02/23/2014)
%%
if length(ytrue)~= length(ypred)
    error('dimension mismatch')
end

ytrue=ytrue(:);
ypred=ypred(:);

n = size(ytrue,1);

%==================================================================================
% check if the labels are all +1 or -1
%==================================================================================
if n~=( sum(ytrue==1)+sum(ytrue==-1) ) || n~=( sum(ypred==1)+sum(ypred==-1) )
    warning(['The entries of ''ytrue'' and ''ypred'' must all be +1 or -1...',...
             'return only the overall accuracy'])
    accuracy  = sum(ytrue==ypred)/n;
    return
end

% overall accuracy
accuracy  = sum(ytrue==ypred)/n;