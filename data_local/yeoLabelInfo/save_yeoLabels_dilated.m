%% save_yeoLabels_dilated.m
% (05/29/2014)
%=========================================================================%
% - save the 'dilated' netowrks (any node that falls within 'radius' mm 
% of a Yeo network is assigned to the network).
% - to be ran after save_yeoLabels.m
%-------------------------------------------------------------------------%
% Outfile: 
%      'yeoLabelInfo_dilated[]mm.mat'
% OutVariables:
%      'roiMNI', 'roiLabel', 'yeoLabels'
%=========================================================================%
%%
clear
purge
%% configure options
fsave=true;

parcellation = 'Grid326'; % {'Grid326','Grid1068','WashU'}

% radius in mm's
radius=5;

outPath=[get_rootdir, '/data_local/yeoLabelInfo/yeo_info_',parcellation, ...
    '_dilated',num2str(radius),'mm.mat'];
outVars={'yeoLabels','roiMNI','roiLabel','mFileName','timeStamp'};
% return
%%
% load output from save_yeoLabels.m
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',parcellation,'.mat'], ...
    'yeoLabels', 'roiMNI')

% path with the YeoPlus label scheme
yeoPlusPath = [get_rootdir,'/data_local/yeoLabelInfo/YeoPlus.hdr'];

% do dilation
roiLabel=mc_NearestNetworkNode(yeoPlusPath, roiMNI, radius);

% round off numerical precision 
roiLabel=round(roiLabel);

mFileName=mfilename('fullpath');
timeStamp=tak_timestamp;
save(outPath,outVars{:})