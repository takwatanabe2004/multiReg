function [w,idx_supp,supp_info,C] = tak_sim_groundTruthWeight_4d_bipartite(nx,ny)
% [w,W,idx_supp] = tak_sim_groundTruthWeight_4d_bipartite(nx,ny)
%=========================================================================%
% - Generate 4d grouth truth weight vector 
%   (4d connectome from 2d node structure)
% - (nx, ny) = # nodes in (x,y)-direction
% - graph consists of a bipartite graph between a pair of node clusters
% - optionally create diffmat for this 4d connectome space
%=========================================================================%
% (07/05/2014)
%%
d = nx*ny; % number of nodes
p = nchoosek(d,2); % number of correlations/edges

%-------------------------------------------------------------------------%
% two clusters of nodes
%-------------------------------------------------------------------------%
% first cluster
idx_nx1 = [6:10];
idx_ny1 = [8:10];
[xx1,yy1]=ndgrid(idx_nx1,idx_ny1);
idx_anom1 = sub2ind([nx,ny],xx1(:),yy1(:))';

% 2nd cluster
idx_nx2 = [3:6];
idx_ny2 = [2:5];
[xx2,yy2]=ndgrid(idx_nx2,idx_ny2);
idx_anom2 = sub2ind([nx,ny],xx2(:),yy2(:))';

supp_info.idx_anom = [idx_anom1, idx_anom2];

%-------------------------------------------------------------------------%
% get indices of node cluster pairs
%-------------------------------------------------------------------------%
[idx_supp,supp_info.mask,supp_info.maskMat]=tak_nodes2edges_2clusters(nx,ny,idx_anom1,idx_anom2);
%% assign ground truth weight on the support
w = zeros(p,1);
w(idx_supp) = 5+5*rand + 1*randn([length(idx_supp),1]);
%%
if nargout==4
    %=====================================================================%
    % create differencing matrix
    %=====================================================================%
    %---------------------------------------------------------------------%
    % lexico-indices of sampled points 
    % (excludes diagonal and upper-triangular coordinates in 4d space)
    %---------------------------------------------------------------------%
    idx_samp = tak_dvec(reshape(1:d^2,[d,d]));

    % 4d adjacency matrix
    adjmat = tak_adjmat_subsampled([nx,ny,nx,ny],idx_samp);
    C = tak_adjmat2incmat(adjmat);
end