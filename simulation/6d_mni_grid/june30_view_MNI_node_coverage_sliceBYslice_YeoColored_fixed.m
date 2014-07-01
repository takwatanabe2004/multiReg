%% june30_view_MNI_node_coverage_sliceBYslice_YeoColored_fixed.m
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
load([get_rootdir,'/data_local/yeoLabelInfo/yeoColorCode.mat'],'yeoColors')
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',GRID, ...
   '_dilated5mm.mat'], 'roiLabel','yeoLabels') %{'roiMNI','roiLabel','yeoLabels'}; 

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
d = size(roiMNI,1); % # of nodes
p = nchoosek(d,2);  % # edges

screenSize=[1 122 1920 855];
%%
figure,set(gcf,'Units','pixels','Position', screenSize)
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

figure,set(gcf,'Units','pixels','Position', screenSize)
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

figure,set(gcf,'Units','pixels','Position', screenSize)
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