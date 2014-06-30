%% june30_view_MNI_node_coverage_sliceBYslice_YeoColored2.m
% (06/30/2014)
%=========================================================================%
% - script builds on june30_view_MNI_node_coverage_sliceBYslice_withIndices.m....
% - plot nodes color-coded with YeoPlus colorcode
%-------------------------------------------------------------------------%
% - view the coverage of the nodes slice-by-slice...from all 3 planes: 
%   coronal, saggital, axial (run the script to see what i mean...)
%=========================================================================%
%%
clear
purge

fsavefig=false;

rootdir = fileparts(mfilename('fullpath'));
GRID='Grid326'; % {'Grid326','Grid1068','WashU'}

% dataPath1=[get_rootdir,'/data_local/designMatrix_FC_',GRID,'.mat'];
% dataVars1={'SubjList'}; % {'X','sex','age'}
% load(dataPath1,dataVars1{:})

dataPath2 = [get_rootdir,'/data_local/yeoLabelInfo/yeo_info_', ...
                GRID,'_dilated5mm.mat'];
dataVars2 = {'roiLabel','yeoLabels'}; %{'roiMNI','roiLabel', 'yeoLabels'}; 
load(dataPath2, dataVars2{:})

%-------------------------------------------------------------------------%
% for some darn reason, there is a node with label '13' in Grid1068...
% set this to 0 (unlabeled)
%-------------------------------------------------------------------------%
if strcmpi(GRID,'Grid1068')
    roiLabel(roiLabel==13)=0;
end

dataPath=[get_rootdir, '/data_local/graphinfo/graph_info_',GRID,'.mat'];
dataVars={'adjmat', 'C', 'coord', 'roiMNI', 'roiMNI_flipx'};
load(dataPath,dataVars{:})
coord.NSIZE
nx=coord.nx;
ny=coord.ny;
nz=coord.nz;
%% plot yeo color scheme
% yeoLabels
% load(['C:\Users\takanori\Documents\MATLAB\BNV_repos\BrainNetViewer',...
%     '\Data\Brain_AAL_Nodes_Edges_edited.mat'])
% yeoColors = [1 1 1; EC.nod.CM(1:12,:)];
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
clear xlen ylen
% savefig([rootdir,'/yeoColorCodes'],'png')
return
%%
d = size(roiMNI,1); % # of nodes
p = nchoosek(d,2);  % # edges
%%
figure,set(gcf,'Units','pixels','Position', [1 122 1920 855])
for x=1:coord.nx
    idx_xslice = find(coord.r(:,1)==x);
    xslice = coord.r(idx_xslice,2:3);
    if strcmpi(GRID,'Grid1068')
        H=subplot(3,4,x);
    else
        H=subplot(2,4,x);
    end
    tak_plot_sim_nodes2d_subplots_yeo_fixed(H,ny,nz,xslice,idx_xslice,...
        roiLabel,yeoColors,0)
end
if fsavefig
    savefig([rootdir,'/MNI_node_coverage_',GRID,'_saggital'],'png')
end

figure,set(gcf,'Units','pixels','Position', [1 122 1920 855])
for y=1:coord.ny
    idx_yslice = find(coord.r(:,2)==y);
    yslice = coord.r(idx_yslice,[1,3]);
    if strcmpi(GRID,'Grid1068')
        H=subplot(3,5,y);
    else
        H=subplot(2,5,y);
    end
    tak_plot_sim_nodes2d_subplots_yeo_fixed(H,nx,nz,yslice,idx_yslice,...
        roiLabel,yeoColors,0)
end
if fsavefig
    savefig([rootdir,'/MNI_node_coverage_',GRID,'_coronal'],'png')
end

figure,set(gcf,'Units','pixels','Position', [1 122 1920 855])
for z=1:coord.nz
    idx_zslice = find(coord.r(:,3)==z);
    zslice = coord.r(idx_zslice,1:2);
    if strcmpi(GRID,'Grid1068')
        H=subplot(2,5,z);
    else
        H=subplot(2,4,z);
    end
    tak_plot_sim_nodes2d_subplots_yeo_fixed(H,nx,ny,zslice,idx_zslice,...
        roiLabel,yeoColors,0)
end
if fsavefig
    savefig([rootdir,'/MNI_node_coverage_',GRID,'_axial'],'png')
end