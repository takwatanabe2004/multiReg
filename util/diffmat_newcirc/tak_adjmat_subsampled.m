function A_subsamp = tak_adjmat_subsampled(ARRAYSIZE, coord)
% A_subsamp = tak_adjmat_subsampled(ARRAYSIZE, coord)
%=========================================================================%
% - create adjacency matrix corresponding to a subsampled version of the 
%   image/volume/connectome.
%-------------------------------------------------------------------------%
% - ARRAYSIZE = [nx,ny,nz...]  # coordinates in each dimension
% - coord = index location of the subsampled coordinates
%         (nsamp x 1) if given in lexicographic index
%         (nsamp x NDIM) if given in (X,Y,Z...) coordinate index
%=========================================================================%
% (06/28/2014)
%%
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

% full adjacency matrix
A = tak_adjmat(ARRAYSIZE);

% subsample adjacency matrix at sampled points
A_subsamp = A(coord,coord);