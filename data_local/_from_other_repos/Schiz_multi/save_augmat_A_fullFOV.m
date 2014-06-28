%% save_augmat_A_fullFOV
% (05/08/2014)
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

mask_nii = load_nii([get_rootdir,'/data_local/brain_mask.nii']);
mask_nii.img=logical(mask_nii.img);

mask_nii=tak_downsample_nii(mask_nii);
mask = mask_nii.img;
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

A_full = A;
timeStamp=tak_timestamp;
mFileName=mfilename;
save([get_rootdir,'/data_local/A_matrix_fullFOV'],'A_full',...
    'nx','ny','nz','timeStamp','mFileName')