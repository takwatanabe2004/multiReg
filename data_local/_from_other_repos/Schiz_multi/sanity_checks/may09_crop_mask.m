%% may09_crop_mask.m
% (05/09/2014)
%=========================================================================%
% - The data augmentation matrix was w.r.t. the entire FOV containing the
%   black background....try again cropiing the background
%=========================================================================%
%%
clear
purge
dir_head = 'C:\Users\takanori\Desktop\Subjects/';
dir_tail = '/session_1/rest_1/coReg/';

load([get_rootdir,'/data_local/A_matrix'],'A')

nii = load_nii(strcat(dir_head,'0040000',dir_tail,'m0wrp1mprage.nii'));

nii.img(isnan(nii.img))=0;
nii.img=double(nii.img);

mask_nii = load_nii([get_rootdir,'/data_local/brain_mask.nii']);
nii=tak_downsample_nii(nii);
mask_nii=tak_downsample_nii(mask_nii);
mask_nii.img=logical(mask_nii.img);

nii_masked=nii;
nii_masked.img = nii.img.*mask_nii.img;

mask = mask_nii.img;

% nii.img=nii.img.*mask;
%%
% tmp1=load_nii(strcat(dir_head,'0040000',dir_tail,'m0wrp1mprage.nii'));
% tmp2=load_nii([get_rootdir,'/data_local/brain_mask.nii']);
% view_nii(tmp1)
% view_nii(tmp2)
% showcs3(tmp1.img)
% showcs3(tmp2.img)
% return
% view_nii(mask_nii)
% view_nii(nii)
%% create (p x 3) matrix indicating the (x,y,z) coordinate of feature vector x
p = sum(mask(:));
[nx,ny,nz]=size(mask);
[xx,yy,zz]=ndgrid(1:nx,1:ny,1:nz);

%=========================================================================%
% - r: indexing of the 3d coordinate
%=========================================================================%
coord.r = zeros(p,3);
coord.r(:,1)=xx(mask);
coord.r(:,2)=yy(mask);
coord.r(:,3)=zz(mask);
coord.nx=nx;
coord.ny=ny;
coord.nz=nz;
coord.xrange = [min(unique(coord.r(:,1))),    max(unique(coord.r(:,1)))];
coord.yrange = [min(unique(coord.r(:,2))),    max(unique(coord.r(:,2)))];
coord.zrange = [min(unique(coord.r(:,3))),    max(unique(coord.r(:,3)))]
%% crop/trim mask
nii_crop = nii;
nii_crop.img = nii_crop.img( coord.xrange(1):coord.xrange(2), ...
                                           coord.yrange(1):coord.yrange(2),...
                                           coord.zrange(1):coord.zrange(2));
mask_nii_crop = mask_nii;
mask_nii_crop.img = mask_nii_crop.img( coord.xrange(1):coord.xrange(2), ...
                                           coord.yrange(1):coord.yrange(2),...
                                           coord.zrange(1):coord.zrange(2));

%=========================================================================%
% adjust nifti coordinate info
%=========================================================================%
mask_nii_crop.hdr.dime.dim=[3,size(mask_nii_crop.img),1,1,1,1];
mask_nii_crop.hdr.hist.originator(1:3)=[...
    mask_nii_crop.hdr.hist.originator(1)-coord.xrange(1)+1,...
    mask_nii_crop.hdr.hist.originator(2)-coord.yrange(1)+1,...
    mask_nii_crop.hdr.hist.originator(3)-coord.zrange(1)+1];

nii_crop.hdr.dime.dim=[3,size(nii_crop.img),1,1,1,1];
nii_crop.hdr.hist.originator(1:3)=[...
    nii_crop.hdr.hist.originator(1)-coord.xrange(1)+1,...
    nii_crop.hdr.hist.originator(2)-coord.yrange(1)+1,...
    nii_crop.hdr.hist.originator(3)-coord.zrange(1)+1];

% return
view_nii(nii_masked);
view_nii(nii_crop);
view_nii(mask_nii);
view_nii(mask_nii_crop);
%% sanity check on the cropped volume
mask_crop=mask( coord.xrange(1):coord.xrange(2), ...
                 coord.yrange(1):coord.yrange(2),...
                 coord.zrange(1):coord.zrange(2));

%=========================================================================%
% same nnz's?
%=========================================================================%
sum(mask(:))
sum(mask_crop(:))

%=========================================================================%
% feature extraction from original nifti and masked nifti returns same result?
%=========================================================================%
x_ver1 = nii.img(mask);
x_ver2 = nii_crop.img(mask_crop);
if isequal(x_ver1,x_ver2)
    disp('success')
else
    error('welp...')
end
% purge
% tak_gui_show_slices(mask)
% tak_gui_show_slices(mask_crop)