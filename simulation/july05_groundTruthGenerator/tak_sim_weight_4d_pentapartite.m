function [w,idx_supp,supp_info,C] = tak_sim_weight_4d_pentapartite(nx,ny)
% [w,W,idx_supp] = tak_sim_weight_4d_quadpartite(nx,ny)
%=========================================================================%
% - Generate 4d grouth truth weight vector 
%   (4d connectome from 2d node structure)
% - (nx, ny) = # nodes in (x,y)-direction
% - graph consists of a 4-partite graph between a pair of node clusters
% - optionally create diffmat for this 4d connectome space
%=========================================================================%
% (07/05/2014)
%%
d = nx*ny; % number of nodes
p = nchoosek(d,2); % number of correlations/edges

%-------------------------------------------------------------------------%
% assign clusters of nodes
%-------------------------------------------------------------------------%
% first cluster
% idx_nx1 = [6:10];
% idx_ny1 = [8:10];
% [xx1,yy1]=ndgrid(idx_nx1,idx_ny1);
% idx_anom1 = sub2ind([nx,ny],xx1(:),yy1(:))';
idx_anom1 = [19,20, 33:36, 48:51, 64,65];

% 2nd cluster
% idx_nx2 = [3:6];
% idx_ny2 = [2:5];
% [xx2,yy2]=ndgrid(idx_nx2,idx_ny2);
% idx_anom2 = sub2ind([nx,ny],xx2(:),yy2(:))';
idx_anom2 = [175:179, 190:194, 205:209];

% 3rd cluster
idx_anom3 = [25:29, 40:44, 55:59];

% 4th cluster
idx_anom4 = [122:126, 108:110, 94, 138:140, 154]+4-ny;

% 5th cluster
idx_anom5 = [168:170,183:185,198:200];

supp_info.idx_anom = [idx_anom1, idx_anom2, idx_anom3, idx_anom4,idx_anom5];

%-------------------------------------------------------------------------%
% get indices of node cluster pairs
%-------------------------------------------------------------------------%
% 1-2
supp_info.idx_list{1,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom1,idx_anom2);
supp_info.idx_list{1,2}='1-2';

% 1-3
supp_info.idx_list{2,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom1,idx_anom3);
supp_info.idx_list{2,2}='1-3';

% 1-4
supp_info.idx_list{3,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom1,idx_anom4);
supp_info.idx_list{3,2}='1-4';

% 1-5
supp_info.idx_list{4,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom1,idx_anom5);
supp_info.idx_list{4,2}='1-5';

% 2-3
supp_info.idx_list{5,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom2,idx_anom3);
supp_info.idx_list{5,2}='2-3';

% 2-4
supp_info.idx_list{6,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom2,idx_anom4);
supp_info.idx_list{6,2}='2-4';

% 2-5
supp_info.idx_list{7,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom2,idx_anom5);
supp_info.idx_list{7,2}='2-5';

% 3-4
supp_info.idx_list{8,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom3,idx_anom4);
supp_info.idx_list{8,2}='3-4';

% 3-5
supp_info.idx_list{9,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom3,idx_anom5);
supp_info.idx_list{9,2}='3-5';

% 4-5
supp_info.idx_list{10,1}=tak_nodes2edges_2clusters(nx,ny,idx_anom4,idx_anom5);
supp_info.idx_list{10,2}='4-5';

%-------------------------------------------------------------------------%
% combine indices 
%-------------------------------------------------------------------------%
idx_supp = [supp_info.idx_list{1,1}; ...
            supp_info.idx_list{2,1}; ...
            supp_info.idx_list{3,1}; ...
            supp_info.idx_list{4,1}; ...
            supp_info.idx_list{5,1}; ...
            supp_info.idx_list{6,1}; ...
            supp_info.idx_list{7,1}; ...
            supp_info.idx_list{8,1}; ...
            supp_info.idx_list{9,1}; ...
            supp_info.idx_list{10,1}];
        
supp_info.mask = false(p,1);
supp_info.mask(idx_supp) = true;
%% assign ground truth weight on the support
w = zeros(p,1);

%-------------------------------------------------------------------------%
% assign weights "cluster-by-cluster"
for i=1:size(supp_info.idx_list,1)
    idx_tmp = supp_info.idx_list{i,1};
    % make weight piecewise constant + perturbation of noise
    w(idx_tmp) = tak_sample_signed_unif([4,10],1) ...
                  + 1*randn([length(idx_tmp),1]);
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