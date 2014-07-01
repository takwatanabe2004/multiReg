%% june30_assign_anomNodeClusters_and_visualize_sliceBYslice.m
% (06/30/2014)
%=========================================================================%
% - script builds on june30_view_MNI_node_coverage_sliceBYslice_YeoColored_fixed.m
% - assign K - clusters of anomalous nodes....these will form the
%   anomalous edge set consisting of a K-partite graphs
% - the ANOM_NODES are highlighted in boxes
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
nx=coord.nx;
ny=coord.ny;
nz=coord.nz;
d = size(roiMNI,1); % # of nodes
p = nchoosek(d,2);  % # edges

screenSize=[1 122 1920 855];
%%
% figure,imexpl,hist(roiLabel, unique(roiLabel)),
% [xx,yy]=hist(roiLabel, unique(roiLabel))

cnt=1;
% assign K anomalous node clusters from each network membership
if strcmp(GRID,'Grid1068')
%     tmp = [9:12, 85:88, 205:206];
%     idx_anomNodesYeoCluster{13} = [tmp, tmp+ny,tmp+2*ny];
elseif strcmp(GRID,'Grid326')
    
    % cluster1: visual (1:BLUE)
%     tmp = [285:286, 291:292, 315:318, 319:322];
%     idx_anomNodeClusters{cnt} = [tmp]; 
%     boxColorList{cnt} = [1 0.2, 1]; % light pink
%     cnt=cnt+1;

    % cluster2: DA (3:GREEN)
    tmp = [311:314, 273:276, 278:281];
    idx_anomNodeClusters{cnt} = [tmp]; 
    boxColorList{cnt} = 'c';
    cnt=cnt+1;
%     
    % cluster3: frotoparietal (6: ORANGE) & default (7: RED)
    tmp = [153:154,213:214,...
           148:151,209:212,269:272,...
           143:144,205:206,265:266,307:310];
    idx_anomNodeClusters{cnt} = [tmp]; 
    boxColorList{cnt} = 'k';
    cnt=cnt+1;
% 
    % cluster4: cerebellum (12:GRAY)
    tmp = [5:9, 40:42, 11:15, 47:49];
    idx_anomNodeClusters{cnt} = [tmp]; 
    boxColorList{cnt} = 'y';
    cnt=cnt+1;
end
Kcluster = length(idx_anomNodeClusters);


% return
%%
figure,set(gcf,'Units','pixels','Position', screenSize)

%=========================================================================%
% (y,z) saggital plane
%=========================================================================%
for x=1:coord.nx
    idx_xslice = find(coord.r(:,1)==x); % <- in lexico index
    xslice = coord.r(idx_xslice,2:3);   % <- in (y,z) coordinate
    %=====================================================================%
    % find if there are anom node cluster members in this slice
    % - keep this in an (K x 1) cell
    %=====================================================================%
    idx_anomNodeSlice=cell(Kcluster,1);
    for icluster=1:Kcluster
        tmp1=idx_anomNodeClusters{icluster};
        tmp2=idx_xslice';
        idx_anomNodeSlice{icluster}=tmp1(ismember(tmp1,tmp2));
    end    
    
    if strcmpi(GRID,'Grid1068')
        H=subplot(3,4,x);
    else
        H=subplot(2,4,x);
    end
    tak_plot_sim_nodes2d_subplots_anomNodeClusters_yeo(H,ny,nz,xslice,idx_xslice,...
        roiLabel,yeoColors,idx_anomNodeSlice,boxColorList,1)
end
% return
if fsavefig
    savefig([rootdir,'/MNI_node_coverage_',GRID,'_saggital'],'png')
end

%=========================================================================%
% (x,z) coronal plane
%=========================================================================%
figure,set(gcf,'Units','pixels','Position', screenSize)
for y=1:coord.ny
    idx_yslice = find(coord.r(:,2)==y);
    yslice = coord.r(idx_yslice,[1,3]);
    %=====================================================================%
    % find if there are anom node cluster members in this slice
    % - keep this in an (K x 1) cell
    %=====================================================================%
    idx_anomNodeSlice=cell(Kcluster,1);
    for icluster=1:Kcluster
        tmp1=idx_anomNodeClusters{icluster};
        tmp2=idx_yslice';
        idx_anomNodeSlice{icluster}=tmp1(ismember(tmp1,tmp2));
    end    
    if strcmpi(GRID,'Grid1068')
        H=subplot(3,5,y);
    else
        H=subplot(2,5,y);
    end
    tak_plot_sim_nodes2d_subplots_anomNodeClusters_yeo(H,nx,nz,yslice,idx_yslice,...
        roiLabel,yeoColors,idx_anomNodeSlice,boxColorList,1)
end
% return
if fsavefig
    savefig([rootdir,'/MNI_node_coverage_',GRID,'_coronal'],'png')
end

%=========================================================================%
% (x,y) axial plane
%=========================================================================%
figure,set(gcf,'Units','pixels','Position', screenSize)
for z=1:coord.nz
    idx_zslice = find(coord.r(:,3)==z);
    zslice = coord.r(idx_zslice,1:2);
    %=====================================================================%
    % find if there are anom node cluster members in this slice
    % - keep this in an (K x 1) cell
    %=====================================================================%
    idx_anomNodeSlice=cell(Kcluster,1);
    for icluster=1:Kcluster
        tmp1=idx_anomNodeClusters{icluster};
        tmp2=idx_zslice';
        idx_anomNodeSlice{icluster}=tmp1(ismember(tmp1,tmp2));
    end   
    if strcmpi(GRID,'Grid1068')
        H=subplot(2,5,z);
    else
        H=subplot(2,4,z);
    end
    tak_plot_sim_nodes2d_subplots_anomNodeClusters_yeo(H,nx,ny,zslice,idx_zslice,...
        roiLabel,yeoColors,idx_anomNodeSlice,boxColorList,1)
end
if fsavefig
    savefig([rootdir,'/MNI_node_coverage_',GRID,'_axial'],'png')
end