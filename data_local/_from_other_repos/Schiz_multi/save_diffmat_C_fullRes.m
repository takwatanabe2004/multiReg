%% save_diffmat_C_fullRes.m
% (05/10/2014)
%=========================================================================%
% - Construct finite differencing matrix "C", as well as adjacency matrix
% - protocode from may08_get_diffmat.m
% - code mostly from t_get_347_graphinfo.m in my old project
%=========================================================================%
%%
clear
purge

mask_nii = load_nii([get_rootdir,'/data_local/brain_mask.nii']);
% mask_nii=tak_downsample_nii(mask_nii); % <- part removed from save_diffmat_C
mask_nii.img=logical(mask_nii.img);
mask = mask_nii.img;
%% extract feature vector via masking
p = sum(mask(:));
[nx,ny,nz]=size(mask);
P = numel(mask);
%% create (p x 3) matrix indicating the (x,y,z) coordinate of x
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
%% apply the same brute force scheme as t_get_347_graphinfo.m in my old project
%% create a nearest-neighbor graph in 3d...do it brute force
adjmat=sparse(p,p);
timeTotal=tic;
tic
for i=1:p
    if mod(i,2000)==0; i/p*100,toc, end;
    
    % find the index-set of the 1st order nearest neighbor
    idx_NN=sum(abs(bsxfun(@minus, coord.r, coord.r(i,:))),2)==1;
    
    % nearest-neighbors
    adjmat(i,idx_NN)=1;
    
    %=====================================================================%
    % sanity check of the nearest neighbor property...
    % - the selected neighboring edge-sets should have l1 distance of 1 
    %   wrt the in (x,y,z) coordinate in question
    %=====================================================================%
    dist_set=abs(sum(bsxfun(@minus,coord.r(idx_NN==1,:),coord.r(i,:)),2));
    if ~all(bsxfun(@eq, dist_set, 1))
        error('argh..')
    end
end
C=tak_adjmat2incmat(adjmat);
%%
mFileName=mfilename;
timeStamp=tak_timestamp;

save([get_rootdir,'/data_local/diffmat_C_fullRes.mat'],'C','adjmat',...
    'mFileName','timeStamp')