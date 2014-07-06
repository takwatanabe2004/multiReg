%% save_6d_MNI_anomNodeClusters.m
% (06/30/2014)
%=========================================================================%
% - code simply from june30_assign_anomNodeClusters_and_visualize_sliceBYslice.m
% - i also added a cell data that contains the "cluster-by-cluster-edge" indices,
%   so that I can assign different "patch-weights" to each cluster-connections
%-------------------------------------------------------------------------%
% - view the coverage of the nodes slice-by-slice...from all 3 planes: 
%   coronal, saggital, axial (run the script to see what i mean...)
%=========================================================================%
%%
clear
purge

GRID='Grid326'; % {'Grid326','Grid1068','WashU'}

rootdir = fileparts(mfilename('fullpath'));

clusterName = 'anomCluster1';

%-------------------------------------------------------------------------%
% save figure?
%-------------------------------------------------------------------------%
fsavefig=true; %  {'true','false'}
figname = [rootdir,'/',clusterName,'_',GRID];
% return

%-------------------------------------------------------------------------%
% save weight vector support?
%-------------------------------------------------------------------------%
fsaveSupport=true; %  {'true','false'} % save weight vector support?
outPath = [rootdir,'/',clusterName,'_',GRID];
outVars = {'wsupp','Wsupp', 'idx_anomEdgeClusters', 'idx_anomNodeClusters', ...
    'timeStamp','mFileName'};
% return
%%
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

if fsavefig
    screenSize=[1 52 1920 1035];
else
    screenSize=[1 111 1920 888];
end
%%
% figure,imexpl,hist(roiLabel, unique(roiLabel)),
% [xx,yy]=hist(roiLabel, unique(roiLabel))

cnt=1;
% assign K anomalous node clusters from each network membership
if strcmp(GRID,'Grid1068')
%     tmp = [9:12, 85:88, 205:206];
%     idx_anomNodesYeoCluster{13} = [tmp, tmp+ny,tmp+2*ny];
elseif strcmp(GRID,'Grid326')
    % cluster: visual (1:violet)
    tmp = [155:158, 159:163,...
           93:96,97:101];
    idx_anomNodeClusters{cnt} = [tmp]; 
    boxColorList{cnt} = [.5 1, 0.2]*1; % light green
    cnt=cnt+1;
    
    % cluster: frontoparietal (2:BLUE)
    tmp = [285:286, 291:292, 315:318, 319:322];
    idx_anomNodeClusters{cnt} = [tmp]; 
    boxColorList{cnt} = [1 0.2, 1]; % light pink
    cnt=cnt+1;

    % cluster: DA (3:GREEN)
    tmp = [311:314, 273:276, 278:281];
    idx_anomNodeClusters{cnt} = [tmp]; 
    boxColorList{cnt} = 'c';
    cnt=cnt+1;
%     
    % cluster: frotoparietal (6: ORANGE) & default (7: RED)
    tmp = [153:154,213:214,...
           148:151,209:212,269:272,...
           143:144,205:206,265:266,308:309];
    idx_anomNodeClusters{cnt} = [tmp]; 
    boxColorList{cnt} = 'k';
    cnt=cnt+1;
% 
    % cluster: cerebellum (12:GRAY)
    tmp = [5:9, 40:42, 11:15, 47:49];
    idx_anomNodeClusters{cnt} = [tmp]; 
    boxColorList{cnt} = 'y';
    cnt=cnt+1;
end
K = length(idx_anomNodeClusters);
%% number of edges in a complete K-partite graph?
for k=1:K
    numNodesList(k) = length(idx_anomNodeClusters{k});
end

nEdges = 0;
for i=1:K-1
    for j=i+1:K
        nEdges = nEdges + numNodesList(i)*numNodesList(j);
    end
end
nEdges
nEdges/p

%% visualize impacted edges in connectome space
% purge
Wsupp = zeros(d,d);
wsupp = [];
idx_anomEdgeClusters = cell(K,1);
cnt=1;
for i=1:K-1
    for j=i+1:K
%         disp('-----------')
%         idx_anomNodeClusters{i}
%         idx_anomNodeClusters{j}
        Wsupp(idx_anomNodeClusters{i},idx_anomNodeClusters{j})=1;
        Wsupp(idx_anomNodeClusters{j},idx_anomNodeClusters{i})=1;
        
        %=================================================================%
        % i also want the indices of the anamolous edges "cluster-by-cluster"
        %=================================================================%
        Wsupp_tmp = zeros(d,d);
        Wsupp_tmp(idx_anomNodeClusters{i},idx_anomNodeClusters{j})=1;
        Wsupp_tmp(idx_anomNodeClusters{j},idx_anomNodeClusters{i})=1;
%         idx_lex6d = tak_nodes2edges_2clusters_new(coord.NSIZE,...
%                           idx_anomNodeClusters{j},idx_anomNodeClusters{i});
        idx_edgeCluster = find(tak_dvec(Wsupp_tmp));
        wsupp = [wsupp; idx_edgeCluster];
        idx_anomEdgeClusters{cnt,1}=idx_edgeCluster;
        cnt=cnt+1;
%         [i,j]
    end
end
imedger(Wsupp)
% wsupp2 = find(tak_dvec(Wsupp));
% isequal(sort(wsupp), sort(wsupp2))
% return
if fsavefig,savefig([figname,'_edgemat'],'png'),end;
% return
%=========================================================================%
% display in sorted coordinate
%=========================================================================%
% circularly shift 1 indices (so "unlabeled" is at the final label index)
roiLabel_circshift=roiLabel-1;
roiLabel_circshift(roiLabel_circshift==-1)=12;
yeoLabels_circshift=circshift(yeoLabels,-1);
[idxsrt,labelCount_circshift] = tak_get_yeo_sort(roiLabel_circshift);

textOption1={'fontweight','b','fontsize',9};
lineOption = {'color','k','linewidth',0.5};

Wsupp_srt = Wsupp(idxsrt,idxsrt);
imedger(Wsupp_srt),axis off
tak_local_linegroups5(gcf,labelCount_circshift,textOption1,yeoLabels_circshift,lineOption)
if fsavefig,savefig([figname,'_edgemat_sorted'],'png'),end;
% return
%% save data
if fsaveSupport
    timeStamp=tak_timestamp;
    mFileName=mfilename;    
    save(outPath,outVars{:})
end
%% display figure
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
    idx_anomNodeSlice=cell(K,1);
    for icluster=1:K
        tmp1=idx_anomNodeClusters{icluster};
        tmp2=idx_xslice';
        idx_anomNodeSlice{icluster}=tmp1(ismember(tmp1,tmp2));
    end    
    
    if strcmpi(GRID,'Grid1068')
        H=subplot(3,4,x);
    else
        H=subplot(3,3,x);
    end
    tak_plot_sim_nodes2d_subplots_anomNodeClusters_yeo(H,ny,nz,xslice,idx_xslice,...
        roiLabel,yeoColors,idx_anomNodeSlice,boxColorList,1)
end
% return
if fsavefig,savefig([figname,'_saggital'],'png'),end;

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
    idx_anomNodeSlice=cell(K,1);
    for icluster=1:K
        tmp1=idx_anomNodeClusters{icluster};
        tmp2=idx_yslice';
        idx_anomNodeSlice{icluster}=tmp1(ismember(tmp1,tmp2));
    end    
    if strcmpi(GRID,'Grid1068')
        H=subplot(3,5,y);
    else
        H=subplot(3,4,y);
    end
    tak_plot_sim_nodes2d_subplots_anomNodeClusters_yeo(H,nx,nz,yslice,idx_yslice,...
        roiLabel,yeoColors,idx_anomNodeSlice,boxColorList,1)
end
% return
if fsavefig,savefig([figname,'_coronal'],'png'),end;

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
    idx_anomNodeSlice=cell(K,1);
    for icluster=1:K
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

if fsavefig,savefig([figname,'_axial'],'png'),end;