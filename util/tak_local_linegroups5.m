function tak_local_linegroups5(h,labelCount,textOption,yeoLabels,lineOption)
% (09/29/2013)
% (03/02/2014) - modified 
% created for tak_exp_plot_result_noreshuffling
%%
nlabels = length(labelCount);
d = sum(labelCount); % # seeds

% draw lines
offset = 1;
offs=0.5; % <- offset for the line needed so the line won't jump on top of the pixels
line( [0 0]+offs, [1 d]+offs, lineOption{:})
line( [0 d]+offs, [0 0]+offs, lineOption{:})
for i = 1:nlabels % -1 since we don't need the line at the bottom right
    ind = offset + labelCount(i) - 1;    
%     line( [ind ind], [1 d], lineOption{:})
%     line( [1 d], [ind ind], lineOption{:})
    line( [ind ind]+offs, [0 d]+offs, lineOption{:})
    line( [0 d]+offs, [ind ind]+offs, lineOption{:})
    offset = ind + 1;
end

% write labels 
offset = 1;
for i = 1:nlabels    
    ind = offset + labelCount(i) - 1;  
    xx = offset+floor((ind-offset)/2);
    text(d+2,xx,yeoLabels{i},textOption{:})
%     text(xx,d+5,yeoLabels{i},textOption{:},'rotation',-35)

    %-------------------------------------------------------------------------%
    % network membership numbers (modified (03/02/2014))
    %-------------------------------------------------------------------------%
    if i==13
        text(-3,xx,'\times',textOption{:},...
            'HorizontalAlignment','right','VerticalAlignment','middle')
        text(xx,d,'\times',textOption{:},...
            'HorizontalAlignment','center','VerticalAlignment','top')
    else
        text(-3,xx,num2str(i),textOption{:},...
            'HorizontalAlignment','right','VerticalAlignment','middle')
        text(xx,d,num2str(i),textOption{:},...
            'HorizontalAlignment','center','VerticalAlignment','top')
    end
    offset = ind + 1;
end

drawnow