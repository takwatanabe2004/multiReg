function A_subsamp = tak_adjmat_subsampled(ARRAYSIZE, coord)
% A_subsamp = tak_adjmat_subsampled(ARRAYSIZE, coord)
%=========================================================================%
% - create adjacency matrix corresponding to a subsampled version of the 
%   image/volume/connectome.
%-------------------------------------------------------------------------%
% - ARRAYSIZE = [nx,ny,nz...]  # coordinates in each dimension
% - coord = lexicographic index location of the subsampled coordinates
%=========================================================================%
% (06/28/2014)
%%
% full adjacency matrix
A = tak_adjmat(ARRAYSIZE);

% subsample adjacency matrix at sampled points
A_subsamp = A(coord,coord);