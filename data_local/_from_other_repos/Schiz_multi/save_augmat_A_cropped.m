%% save_augmat_A_cropped
% (05/09/2014)
%=========================================================================%
% - save sparse matrix A (XYZ x p)
% - p = # features extracted by mask
% - XYZ = # voxels in rectangular FOV in the original volume space of brain
%-------------------------------------------------------------------------%
% reshape(A*x,[X,Y,Z]); % maps feature vector into original 3d brain space
%                       % (think of it as "data augmentation")
% A' * brainVolume(:);  % maps brain volume to feature vector space
%                       % (equivalent to applying mask)
%=========================================================================%
%%
clear
purge

load([get_rootdir,'/data_local/mask_cropped.mat'],'nii_mask_crop','crop_range')
mask = nii_mask_crop.img;
mask_vec=mask(:);

p = sum(mask_vec);
[nx,ny,nz]=size(mask);
A=sparse( numel(mask), p);

cnt=1;
for i=1:size(A,1)
    if mod(i,10e3)==0; i, end;
    if mask_vec(i)        
        A(i,cnt) = 1;
        cnt=cnt+1;
    end
end
isequal(A'*A,speye(p))

timeStamp=tak_timestamp;
mFileName=mfilename;
save([get_rootdir,'/data_local/A_matrix_cropped'],'A',...
    'nx','ny','nz','timeStamp','mFileName')