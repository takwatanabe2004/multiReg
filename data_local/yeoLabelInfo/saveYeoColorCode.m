%% saveYeoColors
% (06/30/2014)
%=========================================================================%
% - save yeo color code scheme
%=========================================================================%
%%
clear all;
purge

fsave=false;

rootdir = fileparts(mfilename('fullpath'));
outPath = [rootdir,'/yeoColorCode.mat'];
outImage= [rootdir,'/yeoColorCode'];
outVars = {'yeoColors','yeoLabels','timeStamp','mFileName'};
%%
% The YeoPlus labels associated with indices 1...13
yeoLabels = {'Unlabeled','Visual', 'Somatomotor', 'Dorsal Attention', ...
    'Ventral Attention', 'Limbic', 'Frontoparietal', 'Default', ...
    'Striatum', 'Amygdala', 'Hippocampus', 'Thalamus', 'Cerebellum'}';
yeoColors = [ ...
    1 1 1; % unlabeled  -  white
    0.8 0 0.9; %  1. visual  -  purple
    0.3 0.3 1; %  2. somatomotor  -  blue
    0 0.7 0; %  3. DA - green
    1 0.85 0.85; %  4. VA - pink (light purple?)
    1 .9 .8; %  5. Limbic - peach/cream
    1 .5 0; %  6. frontoparietal - orange
    1 0 0; %  7. default - red
    0.2 1 1; %  8. striatum - cyan (light blue?)
    1 1 0; %  9. amygdala - yellow
    1 .4 0.4; % 10. hippicampus - light red
    1 0.4 0.7; % 11. thalamus - pink
    .7 .7 .7; % 12. cerebellum - gray
    ];
xlen=5;
ylen=13;
figure,set(gcf,'Units','pixels','Position', [1000 51 800 950])
imagesc(zeros(ylen,xlen))
% colormap([0.85 1 0.85])
colormap([1 1 1])
hold on
axis('off','image');
msize=40;
markerOption={'o', 'MarkerEdgeColor','k','MarkerSize',msize,'linewidth',3};
textOption={'fontsize',msize/2,'fontweight','b','HorizontalAlignment','left',...
            'VerticalAlignment','middle'};
for i=1:length(yeoLabels)
    plot(1,i,markerOption{:},'MarkerFaceColor', yeoColors(i,:))  
    text(xlen/3,i, yeoLabels{i},textOption{:})
    text(1,i, num2str(i-1),textOption{:},'HorizontalAlignment','center')
end
% return
%% save
timeStamp=tak_timestamp;
mFileName=mfilename;

if fsave
    savefig(outImage,'png')
    save(outPath,outVars{:})
end
% return