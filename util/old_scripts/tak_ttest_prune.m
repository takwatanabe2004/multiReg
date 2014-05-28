function [idxSort,pval] = tak_ttest_prune(X, ylabel)
%| [idxSort,pval] = tak_ttest_prune(X, ylabel)
%|------------------------------------------------------------------------------|%
%| Get the paired-sample-t-test ranking of the features between ylabel\in{-1,+1}
%|  Useful for binary classification
%|------------------------------------------------------------------------------|%
%| Input: 
%|  X = (n x p) design matrix
%|  ylabel = {+1, -1} class label.
%|------------------------------------------------------------------------------|%
%| Output: 
%|   idxSort = feature ranked from "most" discrminative to "least" discriminative,
%|             based on the ttest score.
%|------------------------------------------------------------------------------|%
%| 11/25/2012
%| 03/04/2013 -> WARNING!  previous version assumed the feature-matrix was p x n,
%|               here we assume X = (n x p)  
%| 03/30/2013 -> include (sorted) pval as optional output
%%

%|-------- feature selection by filtering  ----------|%
%| Two-sample t-test
% [h,pval] = ttest2(xfeature(:,ylabel==+1)',xfeature(:,ylabel==-1)');
[h,pval] = ttest2(X(ylabel==+1,:),X(ylabel==-1,:));

% Clean out NaNs by setting to 1 (no significance)
pval(isnan(pval))=1;

featurefitness = 1- pval;
[trash, idxSort] = sort(featurefitness,'descend');

pval = pval(idxSort);
