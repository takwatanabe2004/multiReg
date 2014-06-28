function [A,b]=tak_get_augMat_maskMat(ARRAYSIZE,coord)
% tak_function
%=========================================================================%
% - given the subsampled coordinate points "coord" and array dimension 
%   ARRAYSIZE=[nx, ny, ...], create "augmentation matrix" A and 
%   "masking vector/matrix" b for the fft-based Fused Lasso ADMM.
%-------------------------------------------------------------------------%
% - ARRAYSIZE = [nx,ny,nz...]  # coordinates in each dimension
% - coord = index location of the subsampled coordinates
%         (nsamp x 1) if given in lexicographic index
%         (nsamp x NDIM) if given in (X,Y,Z...) coordinate index
%=========================================================================%
% (06/28/2014)
%% parse inputs
% d_subsamp = # subsampled coordinate points
d_subsamp=size(coord,1);

% d_full = # full coordinate points
d_full = prod(ARRAYSIZE);

NDIM = length(ARRAYSIZE);

% if "coord" is given in (x,y,z,...) coordinate index, convert to lexico-index
switch size(coord,2)
    case 2
    coord = sub2ind(ARRAYSIZE, coord(:,1),coord(:,2));
    case 3
    coord = sub2ind(ARRAYSIZE, coord(:,1),coord(:,2),coord(:,3));
    case 4
    coord = sub2ind(ARRAYSIZE, coord(:,1),coord(:,2),coord(:,3), coord(:,4));
    case 6
    coord = sub2ind(ARRAYSIZE, coord(:,1),coord(:,2),coord(:,3), ...
                               coord(:,4),coord(:,5),coord(:,6));
end
%% augmentation matrix
A=sparse(coord, 1:d_subsamp, 1, d_full,d_subsamp);
%% create binary masking matrix to apply the circulant difference matrix
%=========================================================================%
% get Bsupp := support mask of the NDIM-image support
%=========================================================================%
support_mask=reshape(A*true(d_subsamp,1),ARRAYSIZE);

%-------------------------------------------------------------------------%
% create blockdiagonal matrix using "sparse" directly
% (blkdiag is not efficient as it doesn't specialize to "sparse matrix")
%-------------------------------------------------------------------------%
In = eye(NDIM);
D=NDIM*d_full;
bsupp = false(D,1);
for i=1:NDIM
    tmp=circshift(support_mask,In(i,:))-support_mask;
    idx = (1:d_full) + (i-1)*d_full; 
    bsupp(idx) = tmp(:)==0;
end

%-------------------------------------------------------------------------%
% the above support matrix must be composed with the binary circulant mask
% - bsupp: masks the artifacts from the image support
% - bcirc: masks the "wrap-around effect" from the circulant difference matrix
%-------------------------------------------------------------------------%
bcirc=diag(tak_circmask(ARRAYSIZE));
b=bsupp.*bcirc;

%-------------------------------------------------------------------------%
% approach below creates a diagonal matrix...which is memory inefficient...,
% so i went with the above approach of only handling with the diagonal
% components
%-------------------------------------------------------------------------%
% Bsupp = sparse(D,D); % <- initialize empty sparse matrix
% for i=1:NDIM
%     tmp=circshift(support_mask,In(i,:))-support_mask;
%     idx = (1:d_full) + (i-1)*d_full; 
%     Bsupp = Bsupp + sparse(idx(:),idx(:),tmp(:)==0,D,D);
% end
% keyboard
%% old archaic way of obtaining Bsupp....(script kept for my own reference)
%% (suppress cell if you don't care about the history of this script)
% Bx1=circshift(support_mask,[+1  0  0  0  0  0])-support_mask;
% By1=circshift(support_mask,[ 0 +1  0  0  0  0])-support_mask;
% Bz1=circshift(support_mask,[ 0  0 +1  0  0  0])-support_mask;
% Bx2=circshift(support_mask,[ 0  0  0 +1  0  0])-support_mask;
% By2=circshift(support_mask,[ 0  0  0  0 +1  0])-support_mask;
% Bz2=circshift(support_mask,[ 0  0  0  0  0 +1])-support_mask;
% Bx1=tak_spdiag(Bx1(:)==0);
% By1=tak_spdiag(By1(:)==0);
% Bz1=tak_spdiag(Bz1(:)==0);
% Bx2=tak_spdiag(Bx2(:)==0);
% By2=tak_spdiag(By2(:)==0);
% Bz2=tak_spdiag(Bz2(:)==0);
% 
% % blkdiag can be slow for large sparse matrices
% Bsupp=...
% [          Bx1, sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full); ...
%  sparse(d_full,d_full),           By1, sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full); ...
%  sparse(d_full,d_full), sparse(d_full,d_full),           Bz1, sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full); ...
%  sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full),           Bx2, sparse(d_full,d_full), sparse(d_full,d_full); ...
%  sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full),           By2, sparse(d_full,d_full); ...
%  sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full), sparse(d_full,d_full),           Bz2];
% 
%-------------------------------------------------------------------------%
% the above support matrix must be composed with the binary circulant mask
% - Bsupp: masks the artifacts from the image support
% - Bcirc: masks the "wrap-around effect" from the circulant difference matrix
%-------------------------------------------------------------------------%
% Bcirc=tak_circmask(ARRAYSIZE);
% B=Bsupp*Bcirc;
% b=logical(full(diag(B)));
% keyboard