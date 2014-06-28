%% save_featureMatrix_sMRI_dsamp4.m
% (05/08/2014)
%=========================================================================%
% save design matrix at downsampled resolution (dsamp by factor of 4)
%=========================================================================%
%%
clear
purge
dir_head = 'C:\Users\takanori\Desktop\Subjects/';
dir_tail = '/session_1/rest_1/coReg/';

SubjDir = S_subjDir_Schiz_COBRE_censor;
SubjList = SubjDir(:,1);
n = length(SubjList);


%%
% nii1 = load_nii(strcat(dir_head,'0040000',dir_tail,'m0wrp1mprage.nii'));
mask = load_nii([get_rootdir,'/data_local/brain_mask.nii']);
% return
% view_nii(nii1)
% view_nii(mask)

% nii2=tak_downsample_nii(nii1);
mask=tak_downsample_nii(tak_downsample_nii(mask));
% return
% view_nii(nii2)
% view_nii(mask)
% sum((logical(mask.img)))
mask =logical(mask.img);
p=sum(mask(:));
%%
X = zeros(n,p);
for idx=1:n
    idx
    nii = load_nii(strcat(dir_head,SubjList{idx},dir_tail,'m0wrp1mprage.nii'));
    nii=tak_downsample_nii(tak_downsample_nii(nii));
    x = nii.img(mask);
%     keyboard
    X(idx,:) = x(:)';
end
%%
y = cell2mat(SubjDir(:,2))
timeStamp=tak_timestamp;
mFileName=mfilename;
save([get_rootdir,'/data_local/sMRI_design_censor_dsamp4.mat'],'X','y','mFileName','timeStamp')
%%
% purge
% Xstd=zscore(X);
% figure,imagesc(X)
% figure,imagesc(Xstd)