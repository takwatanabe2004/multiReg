%% may06_test_downsample_nii.m
% (05/06/2014)
%=========================================================================%
% try to downsample nii volume
%=========================================================================%
%%
clear
purge
dir_head = 'C:\Users\takanori\Desktop\Subjects/';
dir_tail = '/session_1/rest_1/coReg/';

nii1 = load_nii(strcat(dir_head,'0040000',dir_tail,'m0wrp1mprage.nii'));
mask = load_nii([get_rootdir,'/data_local/brain_mask.nii']);
view_nii(nii1)
view_nii(mask)

nii2=tak_downsample_nii(nii1);
mask=tak_downsample_nii(mask);
view_nii(nii2)
view_nii(mask)
%%
sum(nii2.img(logical(mask.img)))
x = nii2.img(logical(mask.img));