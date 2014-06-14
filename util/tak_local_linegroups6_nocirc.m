function tak_local_linegroups6_nocirc(h,labelCount,yeoLabels,textOption,lineOption)
% (06/11/2014)
% created rom tak_local_linegroups6, but without circshifting the yeolabels
%%
%-------------------------------------------------------------------------%
% added default settings
%-------------------------------------------------------------------------%
if nargin==3
    textOption={'fontweight','b','fontsize',12};
    lineOption = {'color','k','linewidth',0.5};
end
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
    if i==1
        text(-3,xx,'\times',textOption{:},...
            'HorizontalAlignment','right','VerticalAlignment','middle')
        text(xx,d,'\times',textOption{:},...
            'HorizontalAlignment','center','VerticalAlignment','top')
    else
        text(-3,xx,num2str(i-1),textOption{:},...
            'HorizontalAlignment','right','VerticalAlignment','middle')
        text(xx,d,num2str(i-1),textOption{:},...
            'HorizontalAlignment','center','VerticalAlignment','top')
    end
    offset = ind + 1;
end

drawnow