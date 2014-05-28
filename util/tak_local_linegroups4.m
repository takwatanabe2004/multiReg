function tak_local_linegroups4(h,labelCount,textOption,yeoLabels,lineOption)
% (09/29/2013)
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
%     text(d+2,xx,yeoLabels{i},textOption{:})
%     text(xx,d+5,yeoLabels{i},textOption{:},'rotation',-35)
    FU=9;
    fhorz=12.5;
    fvert=11.5;
    if i==8
        text(-1,xx,num2str(i),'fontweight','b','fontsize',fhorz,...
            'HorizontalAlignment','right','VerticalAlignment','middle')
        text(xx,d,num2str(i),'fontweight','b','fontsize',fvert,...
            'HorizontalAlignment','center','VerticalAlignment','top')
    elseif i==9
        text(-15,xx,num2str(i),'fontweight','b','fontsize',fhorz,...
            'HorizontalAlignment','right','VerticalAlignment','middle')
        text(xx-2,d+FU,num2str(i),'fontweight','b','fontsize',fvert,...
            'HorizontalAlignment','center','VerticalAlignment','top')
    elseif i==10
        text(-1,xx,num2str(i),'fontweight','b','fontsize',fhorz,...
            'HorizontalAlignment','right','VerticalAlignment','middle')
        text(xx,d,num2str(i),'fontweight','b','fontsize',fvert,...
            'HorizontalAlignment','center','VerticalAlignment','top')
    elseif i==11
        text(-15,xx,num2str(i),'fontweight','b','fontsize',fhorz,...
            'HorizontalAlignment','right','VerticalAlignment','middle')
        text(xx+2,d+FU,num2str(i),'fontweight','b','fontsize',fvert,...
            'HorizontalAlignment','center','VerticalAlignment','top')
%     elseif i==12
%         text(-13,xx+2,num2str(i),'fontweight','b','fontsize',fhorz,...
%             'HorizontalAlignment','right','VerticalAlignment','middle')
%         text(xx+3,d+FU,num2str(i),'fontweight','b','fontsize',fvert,...
%             'HorizontalAlignment','center','VerticalAlignment','top')
    elseif i==13
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
%     if i==11
%         text(d+2,xx,num2str(i-1),textOption{:})
%         text(xx-14,d+10,num2str(i-1),textOption{:})
%     elseif i==12
%         text(d+2,xx,num2str(i-1),textOption{:})
%         text(xx-6,d+10,num2str(i-1),textOption{:})
%     else
%         text(d+2,xx,num2str(i-1),textOption{:})
%         text(xx-4,d+10,num2str(i-1),textOption{:})
%     end
    offset = ind + 1;
end

drawnow