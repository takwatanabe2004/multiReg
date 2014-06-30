%% june29_plot_timeSeries_roiTC
% (06/29/2014)
%=========================================================================%
% - Plot actual timeseries from subject from NKI
% - visually study the spatial & temporal correlation pattern
%   (notice the periodic structure in the correlation matrix across time)
%=========================================================================%
%%
clear
purge

GRID = 'Grid1068'; % {'Grid326','Grid1068','WashU'}
% nodeIdx = 115:117;

idx_start=305;
len=2;
nodeIdx = idx_start:idx_start+len;

dirHead=[get_rootdir,'/data_local/0100451/'];

%% concatenated
% load([dirHead,'/Concat_645_1400/s6_',GRID,'_p20f0b_p35mask\s6_',GRID,'_p20f0b_p35mask_roiTC.mat'])
% figure,imexpb,tplot(roiTC(:,nodeIdx))
%% low res, small # samples
% load([dirHead,'/RfMRI_mx_1400/s6_',GRID,'_p20f0b_p35mask/s6_',GRID,'_p20f0b_p35mask_roiTC.mat'])
% figure,imexpb,tplot(roiTC(:,nodeIdx))
%     imcovvl(corr(roiTC'))
%     imcovvl(corr(roiTC))
%% high res, large # samples
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',GRID,'.mat'],'roiMNI')
load([dirHead,'/RfMRI_mx_645/s6_',GRID,'_p20f0b_p35mask/s6_',GRID,'_p20f0b_p35mask_roiTC.mat'])
figure,imexpb,tplot(roiTC(:,nodeIdx)),
%     imcovvr(corr(roiTC(:,nodeIdx))),%imcovvr(corr(roiTC(:,nodeIdx)'))
    imcovvl(corr(roiTC'))
%     imcovvl(corr(roiTC))
    imcovvl(corr(roiTC)),caxis(.5*[-1 1])
%% box plot
coord=tak_discretize_coord(roiMNI);
zlist = unique(coord(:,3));
nz = max(zlist);
for z=1:nz
    zLen(z) = sum( coord(:,3)==z);
end
zOffSets=[0 cumsum(zLen)];
offset=0;
imcovvl(corr(roiTC)),caxis(.5*[-1 1])
for z=1:z
    offset = zOffSets(z);
    len = zLen(z);
    tak_plot_box(offset+[1,len,1,len],'k',3)
end
% t
