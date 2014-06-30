%% june30_view_MNI_node_coverage_sliceBYslice_withIndices.m
% (06/30/2014)
%=========================================================================%
% - script builds on june29_view_MNI_node_coverage_sliceBYslice.m....
%-------------------------------------------------------------------------%
% - view the coverage of the nodes slice-by-slice...from all 3 planes: 
%   coronal, saggital, axial (run the script to see what i mean...)
%=========================================================================%
%%
clear
purge

rootdir = fileparts(mfilename('fullpath'));
GRID='Grid326'; % {'Grid326','Grid1068','WashU'}

% dataPath1=[get_rootdir,'/data_local/designMatrix_FC_',GRID,'.mat'];
% dataVars1={'SubjList'}; % {'X','sex','age'}
% load(dataPath1,dataVars1{:})

% dataPath2 = [get_rootdir,'/data_local/yeoLabelInfo/yeo_info_', ...
%                 GRID,'_dilated5mm.mat'];
% dataVars2 = {'roiMNI'}; %{'roiLabel', 'yeoLabels'}; 
% load(dataPath2, dataVars2{:})

dataPath=[get_rootdir, '/data_local/graphinfo/graph_info_',GRID,'.mat'];
dataVars={'adjmat', 'C', 'coord', 'roiMNI', 'roiMNI_flipx'};
load(dataPath,dataVars{:})
coord.NSIZE
%%
d = size(roiMNI,1); % # of nodes
p = nchoosek(d,2);  % # edges
%%
figure,set(gcf,'Units','pixels','Position', [1 122 1920 855])
cnt=1;
for x=1:coord.nx
    idx_xslice = (coord.r(:,1)==x);
    xslice = coord.r(idx_xslice,2:3);
    if strcmpi(GRID,'Grid1068')
        H=subplot(3,4,x);
    else
        H=subplot(2,4,x);
    end
    tak_plot_sim_nodes2d_subplots3(H,coord.ny,coord.nz,xslice,cnt)
    cnt=cnt+length(xslice);
%     return
end

figure,set(gcf,'Units','pixels','Position', [1 122 1920 855])
cnt=1;
for y=1:coord.ny
    idx_yslice = (coord.r(:,2)==y);
    yslice = coord.r(idx_yslice,[1,3]);
    if strcmpi(GRID,'Grid1068')
        H=subplot(3,5,y);
    else
        H=subplot(2,5,y);
    end
    tak_plot_sim_nodes2d_subplots3(H,coord.nx,coord.nz,yslice,cnt)
    cnt=cnt+length(yslice);
end

figure,set(gcf,'Units','pixels','Position', [1 122 1920 855])
cnt=1;
for z=1:coord.nz
    idx_zslice = (coord.r(:,3)==z);
    zslice = coord.r(idx_zslice,1:2);
    if strcmpi(GRID,'Grid1068')
        H=subplot(2,5,z);
    else
        H=subplot(2,4,z);
    end
    tak_plot_sim_nodes2d_subplots3(H,coord.nx,coord.ny,zslice,cnt)
    cnt=cnt+length(zslice);
end
