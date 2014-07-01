function imconnYeoTriagr(w,GRID,flagLineText)
% imconnYeoTriagr(w,GRID,flagLineText)
%=========================================================================%
% - reshape vectorized connectome into symmetric matrix, and reorder
%   weight vector according to yeoNetworkScheme 
%   (thus, code assumes we're dealing in the real MNI space with a real
%    grid parcellation scheme)
% - code is "heavy-handed", and i don't intend to make this a flexible
%   code...just a convenience script for myself
%-------------------------------------------------------------------------%
% w: (p x 1) vector to be displayed as (d x d) symmetric matrix
% GRID: {'Grid326', 'Grid1068', 'Washu'}
% flagLineText: display text display network name on the side (default: off)
%=========================================================================%
% (07/01/2014)
%%
figure,imexpr
if ~exist('GRID','var') || isempty(GRID)
    GRID = 'Grid326'; % <- I mostly focus on this...
end

if ~exist('flagLineText','var') || isempty(flagLineText)
    flagLineText = 0;
end

textOption={'fontweight','b','fontsize',11};
lineOption = {'color','k','linewidth',1};
titleOption = {'Interpreter','none','fontweight','b','fontsize',16};
%% parse yeo info
%-------------------------------------------------------------------------%
% load data (this is really heavy-handed...but whatever....
%-------------------------------------------------------------------------%
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',GRID, ...
   '_dilated5mm.mat'], 'roiLabel','yeoLabels') %{'roiMNI','roiLabel','yeoLabels'}; 

%-------------------------------------------------------------------------%
% for some darn reason, there is a node with label '13' in Grid1068...
% set this to 0 (unlabeled)
%-------------------------------------------------------------------------%
if strcmpi(GRID,'Grid1068')
    roiLabel(roiLabel==13)=0;
end

% circularly shift 1 indices (so "unlabeled" is at the final label index)
roiLabel=roiLabel-1;
roiLabel(roiLabel==-1)=12;
yeoLabels=circshift(yeoLabels,-1);
[idxsrt,labelCount] = tak_get_yeo_sort(roiLabel);
%%
wmat = tak_dvecinv(w,0);
wmat=wmat(idxsrt,idxsrt);
imtriag(wmat), axis off,%colorbar off
tak_line_yeoGroups(labelCount,lineOption)
tak_text_yeoGroups(labelCount,textOption,yeoLabels,flagLineText)

[nnz_lower, sp_level] = tak_nnz_lower(wmat);
titleStr = sprintf('%g nonzeroes (%g%% sparsity level)', nnz_lower,sp_level*100);
title([inputname(1), ' --- ', titleStr],titleOption{:})
drawnow
end %% end of main function

%% private function
function tak_line_yeoGroups(labelCount,lineOption)
%=========================================================================%
% tooks this from tak_local_linegroups5.m
%=========================================================================%
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
end
%%
function tak_text_yeoGroups(labelCount,textOption,yeoLabels,flagLineText)
% write labels 
nlabels = length(labelCount);
d = sum(labelCount); % # seeds
offset = 1;
for i = 1:nlabels    
    ind = offset + labelCount(i) - 1;  
    xx = offset+floor((ind-offset)/2);
    if flagLineText
        text(d+2,xx,yeoLabels{i},textOption{:})
    end
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
end