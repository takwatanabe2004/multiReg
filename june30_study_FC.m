%% mfileName
% (06/30/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge

GRID = 'Grid326'; % {'Grid326','Grid1068','WashU'}

load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',GRID,'.mat'],'roiMNI')
load([get_rootdir,'/data_local/designMatrix_FC_',GRID,'.mat'],'X')

imconnl(X(1,:))
Xmean = mean(X)';
imconnl(Xmean,1),caxis(caxis/2)
%% box plot z-clusters
coord=tak_discretize_coord(roiMNI);
zlist = unique(coord(:,3));
nz = length(zlist);
for z=1:nz
    zLen(z) = sum( coord(:,3)==z);
end
zoffsetList=[0 cumsum(zLen)];

imconnl(Xmean,1),caxis(caxis/2)
for z=1:z
    offset = zoffsetList(z);
    len = zLen(z);    
%     tak_plot_box(offset+[1,len,1,len],'k',2)
    OFFSET = [-0.5,0.5,-0.5,0.5]; % <- included so that box-corners touches each other
    tak_plot_box(offset+OFFSET+[1,len,1,len],'k',2)
end
%% hierarchical box plot over y and z-clusters
% purge
imconnl(Xmean,1),caxis(caxis/2)
for z=1:z
    offsetz = zoffsetList(z);
    len = zLen(z);
%     tak_plot_box(offsetz+[1,len,1,len],'k',2)
    OFFSETZ = [-0.5,0.5,-0.5,0.5]; % <- included so that box-corners touches each other
    tak_plot_box(offsetz+OFFSETZ+[1,len,1,len],'k',2)
    
    %=====================================================================%
    % inner y-cluster
    %=====================================================================%
    zIDX = coord(:,3)==z;
    ylist=unique( coord(zIDX,2) );
%     yLen = length(ylist);
    for y=1:length(ylist)
        yLen(y) = sum(coord(zIDX,2)==ylist(y));
    end
    yoffsetList=[0 cumsum(yLen)];
    
    for y=1:length(ylist)
        offsety = yoffsetList(y)+offsetz;
        len = yLen(y);
%         tak_plot_box(offsety+[1,len,1,len],'m',2)
        OFFSETY = [-0.5,0.5,-0.5,0.5]; % <- included so that box-corners touches each other
        tak_plot_box(offsety+OFFSETY+[1,len,1,len],'m',2)
    end
    clear yLen
end
zoom on