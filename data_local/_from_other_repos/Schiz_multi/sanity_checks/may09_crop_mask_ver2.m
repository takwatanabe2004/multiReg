%% may09_crop_mask_ver2.m
% (05/09/2014)
%=========================================================================%
% - The data augmentation matrix was w.r.t. the entire FOV containing the
%   black background....try again cropiing the background
%-------------------------------------------------------------------------%
% - this "version" uses the mask saved from save_brain_mask_cropped.m
%=========================================================================%
%%
clear
purge
dir_head = 'C:\Users\takanori\Desktop\Subjects/';
dir_tail = '/session_1/rest_1/coReg/';

load([get_rootdir,'/data_local/mask_cropped.mat'],'nii_mask_crop','crop_range')
nii = load_nii(strcat(dir_head,'0040000',dir_tail,'m0wrp1mprage.nii'));

nii.img(isnan(nii.img))=0;
nii.img=double(nii.img);

nii_mask = load_nii([get_rootdir,'/data_local/brain_mask.nii']);
nii=tak_downsample_nii(nii);
nii_mask=tak_downsample_nii(nii_mask);
nii_mask.img=logical(nii_mask.img);

nii_masked=nii;
nii_masked.img = nii.img.*nii_mask.img;
%% crop/trim mask
nii_crop = nii_mask_crop;
nii_crop.img = nii.img(crop_range.x(1):crop_range.x(2), ...
                       crop_range.y(1):crop_range.y(2),...
                       crop_range.z(1):crop_range.z(2));

% return
view_nii(nii_masked);
view_nii(nii_crop);
view_nii(nii_mask);
view_nii(nii_mask_crop);
%% sanity check on the cropped volume

%=========================================================================%
% same nnz's?
%=========================================================================%
sum(nii_mask.img(:))
sum(nii_mask_crop.img(:))

%=========================================================================%
% feature extraction from original nifti and masked nifti returns same result?
%=========================================================================%
x_ver1 = nii.img(nii_mask.img);
x_ver2 = nii_crop.img(nii_mask_crop.img);
if isequal(x_ver1,x_ver2)
    disp('success')
else
    error('welp...')
end
% purge
% tak_gui_show_slices(mask)
% tak_gui_show_slices(mask_crop)