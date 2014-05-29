%% save_yeoLabels.m
% (05/29/2014)
%=========================================================================%
% Save MIN coordinate information and the corresponding yeoLabels
%-------------------------------------------------------------------------%
% Outfile: 
%      'yeoLabelInfo.mat'
% OutVariables:
%      'roiMNI', 'roiLabel', 'yeoLabels'
%=========================================================================%
%%
clear
purge
fsave=true;
%% 
parcellation = 'WashU'; % {'Grid326','Grid1068','WashU'}
outVars= {'roiLabel', 'yeoLabels','roiMNI','mFileName','timeStamp'};
outPath=[get_rootdir, '/data_local/yeoLabelInfo/yeo_info_',parcellation,'.mat']
% return
%% load relevant data
%=========================================================================%
% Get node coordinates in MNI space
%=========================================================================%
nodeInfoPath = [get_rootdir,'/data_local/node_info/s6_',parcellation,...
    '_p20f0b_p35mask_parameters.mat']
load(nodeInfoPath,'parameters');

% the coordinates of the roi's in MNI space (the centroid coordinate)
roiMNI = parameters.rois.mni.coordinates;
% return

%=========================================================================%
% Get Yeo Info
%=========================================================================%
% path with the YeoPlus label scheme
yeoPlusPath = [get_rootdir,'/data_local/yeoLabelInfo/YeoPlus.hdr']

%-------------------------------------------------------------------------%
% The labels associated with each roi's (YeoPlus network).
% "round" is necessary since there is some roundoff error in the label file
%-------------------------------------------------------------------------%
roiLabel = mc_network_lookup(yeoPlusPath, roiMNI);
roiLabel = round(roiLabel(:,4));    

% The YeoPlus labels associated with indices 1...13
yeoLabels = {'Unlabeled','Visual', 'Somatomotor', 'Dorsal Attention', ...
    'Ventral Attention', 'Limbic', 'Frontoparietal', 'Default', ...
    'Striatum', 'Amygdala', 'Hippocampus', 'Thalamus', 'Cerebellum'}';
%% save 
mFileName=mfilename;
timeStamp=tak_timestamp;
if fsave
    save(outPath,outVars{:})
end