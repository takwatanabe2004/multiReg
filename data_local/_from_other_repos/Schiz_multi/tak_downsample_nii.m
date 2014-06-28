function nii_out = tak_downsample_nii(nii)
% nii_out = tak_downsample_nii(nii)
% (05/06/2014)
%=========================================================================%
% - ultra cheap way to downsample a nii struct by factor of 2
% - the "true" nifti convention most likley not followed properly...
%   i think i can get away with this, but fix if necessary in the future.
%=========================================================================%
%%
nii_out=nii;
nii_out.img = nii_out.img(1:2:end,1:2:end,1:2:end);
% nii2.hdr.dime.pixdim= [1 3 3 3 0 0 0 0];
nii_out.hdr.dime.pixdim= [1 nii.hdr.dime.pixdim(2:4)*2 0 0 0 0];
% nii2.hdr.dime.dim=[3,61,73,61,1,1,1,1];
nii_out.hdr.dime.dim=[3,size(nii_out.img),1,1,1,1];
nii_out.hdr.hist.originator=round(nii.hdr.hist.originator/2);