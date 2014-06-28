%% save_brain_mask_cropped_fullRes.m
% (05/10/2014)
%=========================================================================%
% - The data augmentation matrix was w.r.t. the entire FOV containing the
%   black background....try again cropiing the background
%-------------------------------------------------------------------------%
% - code from may09_crop_mask.m
%=========================================================================%
%%
clear
purge
dir_head = 'C:\Users\takanori\Desktop\Subjects/';
dir_tail = '/session_1/rest_1/coReg/';

mask = load_nii([get_rootdir,'/data_local/brain_mask.nii']);
% mask=tak_downsample_nii(mask);
mask.img=logical(mask.img);
%% create (p x 3) matrix indicating the (x,y,z) coordinate of feature vector x
p = sum(mask.img(:));
[nx,ny,nz]=size(mask.img);
[xx,yy,zz]=ndgrid(1:nx,1:ny,1:nz);

%=========================================================================%
% - r: indexing of the 3d coordinate
%=========================================================================%
coord.r = zeros(p,3);
coord.r(:,1)=xx(mask.img);
coord.r(:,2)=yy(mask.img);
coord.r(:,3)=zz(mask.img);
coord.nx=nx;
coord.ny=ny;
coord.nz=nz;

%=========================================================================%
% range of coordinates to include within the crop
%=========================================================================%
crop_range.x= [min(unique(coord.r(:,1))),    max(unique(coord.r(:,1)))];
crop_range.y= [min(unique(coord.r(:,2))),    max(unique(coord.r(:,2)))];
crop_range.z= [min(unique(coord.r(:,3))),    max(unique(coord.r(:,3)))];
%% crop/trim mask
mask_crop = mask;
mask_crop.img = mask_crop.img( crop_range.x(1):crop_range.x(2), ...
                               crop_range.y(1):crop_range.y(2),...
                               crop_range.z(1):crop_range.z(2));

%=========================================================================%
% adjust nifti coordinate info
%=========================================================================%
mask_crop.hdr.dime.dim=[3,size(mask_crop.img),1,1,1,1];
mask_crop.hdr.hist.originator(1:3)=[...
    mask_crop.hdr.hist.originator(1)-crop_range.x(1)+1,...
    mask_crop.hdr.hist.originator(2)-crop_range.y(1)+1,...
    mask_crop.hdr.hist.originator(3)-crop_range.z(1)+1];
view_nii(mask);
view_nii(mask_crop);

nii_mask_crop=mask_crop;

timeStamp=tak_timestamp;
mFileName=mfilename;
% save([get_rootdir,'/data_local/mask_cropped_fullRes.mat'],'nii_mask_crop','crop_range',...
%     'timeStamp','mFileName')
