function [idx_sort,labelCount,idx_sortvec] = tak_get_yeo_sort(roiLabel)
%|------------------------------------------------------------------------------|%
%|------------------------------------------------------------------------------|%
%| Example:
%   [idx_sort,labelCount] = tak_get_yeo_sort(roiLabel)
%   imedge(E_HC)
%   tak_box_covGroups(gcf, labelCount, yeoLabels) 
%%
[trash,idx_sort] = sort(roiLabel);        
labelCount = hist(roiLabel, [min(roiLabel):max(roiLabel)]);

% get the sorting vector for the flattened version of matrix
d=length(roiLabel);
p=d*(d-1)/2;

idx_sortvec = tak_dvecinv(1:p,0);
idx_sortvec = reshape(idx_sortvec,[d,d]);
idx_sortvec = idx_sortvec(idx_sort,idx_sort);
idx_sortvec = tak_dvec(idx_sortvec);
