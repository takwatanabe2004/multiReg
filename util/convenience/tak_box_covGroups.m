function tak_box_covGroups(h, labelCount, Labels, flagText)
% Example:
% load yeo_info347.mat yeoLabels roiMNI roiLabel
% [idx_sort,labelCount] = tak_get_yeo_sort(roiLabel);
% imedge(v2matsort)
% tak_box_covGroups(gcf, labelCount, yeoLabels)
%%
if(~exist('flagText','var')||isempty(flagText)), flagText = false; end

% set(h);
hold on

nlabels = length(labelCount);

ind1 = 1;
for i = 1:nlabels    
    ind2 = ind1 + labelCount(i) - 1;    
    tak_plot_box([ind1,ind2,ind1,ind2])
    if flagText
         text(ind1,ind1,Labels{i},'fontweight','b','fontsize',16)
    end
    ind1 = ind2 + 1;
end
drawnow