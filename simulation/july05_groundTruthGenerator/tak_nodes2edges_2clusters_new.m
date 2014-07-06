function [idx,maskVec,maskMat] = tak_nodes2edges_2clusters_new(NSIZE,idx_cluster1,idx_cluster2)
% [idx,maskVec,maskMat] = tak_nodes2edges_2clusters_new(NSIZE,idx_cluster1,idx_cluster2)
%=========================================================================%
% - same as tak_nodes2edges_2clusters.m, but nodes are in NSIZE space
%   (in old script, nodes are restricted to 2d space...should fix this in future)
%-------------------------------------------------------------------------%
% - map indices of nodes to edge indices, where we assume that we're given
%   the node locations in lexicographic-coordinates of 2 clusters...
% - we assume this 2 clusters form a complete bipartite graph
%=========================================================================%
% (07/05/2014)
%%
d=prod(NSIZE);
p=nchoosek(d,2);
maskMat=false(d,d);
for ii=1:length(idx_cluster1)
    idx1 = idx_cluster1(ii);
    for jj=1:length(idx_cluster2)
        idx2=idx_cluster2(jj);
        maskMat(idx1,idx2)=true;
        maskMat(idx2,idx1)=true; %<- for symmetry
    end
end

% lower triangular mask
maskLower=tril(true(d),-1);
% keyboard
maskMat=maskMat&maskLower;

% tmp=reshape(1:d^2,[d,d]);
tmp=tak_dvecinv(1:p,0);
idx=tmp(maskMat);

maskVec=false(p,1);
maskVec(idx)=true;
% keyboard