function [w,idx_supp,supp_info,C] = tak_sim_weight_4d_tripartite(nx,ny)
% [w,idx_supp,supp_info,C] = tak_sim_weight_4d_tripartite(nx,ny)
%=========================================================================%
% - Generate 4d grouth truth weight vector 
%   (4d connectome from 2d node structure)
% - (nx, ny) = # nodes in (x,y)-direction
% - graph consists of a K-partite graph between a pair of node clusters
% - optionally create diffmat for this 4d connectome space
%=========================================================================%
% (07/05/2014)
%%
if nargin==0
    nx=11;
    ny=11;
end
d = nx*ny; % number of nodes
p = nchoosek(d,2); % number of correlations/edges

%-------------------------------------------------------------------------%
% assign clusters of nodes
%-------------------------------------------------------------------------%
% first cluster
idx_anom1 = [13:15, 24:26, 35:37, 46:48];

% 2nd cluster
idx_anom2 = [92:97, 103:108];

% 3rd cluster
idx_anom3 = [40:44, 30:32, 20, 52:54,64];

supp_info.idx_anom = [idx_anom1, idx_anom2, idx_anom3];

%-------------------------------------------------------------------------%
% get indices of node cluster pairs
%-------------------------------------------------------------------------%
% 1-2
supp_info.idx_list{1,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom1,idx_anom2);
supp_info.idx_list{1,2}='1-2';

% 1-3
supp_info.idx_list{2,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom1,idx_anom3);
supp_info.idx_list{2,2}='1-3';

% 2-3
supp_info.idx_list{3,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom2,idx_anom3);
supp_info.idx_list{3,2}='2-3';

% combine indices 
idx_supp = [supp_info.idx_list{1,1}; ...
            supp_info.idx_list{2,1}; ...
            supp_info.idx_list{3,1}];
        
supp_info.mask = false(p,1);
supp_info.mask(idx_supp) = true;
supp_info.nx=nx;
supp_info.ny=ny;
%% assign ground truth weight on the support
w = zeros(p,1);
%-------------------------------------------------------------------------%
% assign weights "cluster-by-cluster"
weightList = [8, 4, -10];

for i=1:size(supp_info.idx_list,1)
    idx_tmp = supp_info.idx_list{i,1};
    % make weight piecewise constant + perturbation of noise
    w(idx_tmp) = weightList(i) ...
                  + 0.5*randn([length(idx_tmp),1]);
%     w(idx_tmp) = tak_sample_signed_unif([4,10],1) ...
%                   + 1*randn([length(idx_tmp),1]);
end
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