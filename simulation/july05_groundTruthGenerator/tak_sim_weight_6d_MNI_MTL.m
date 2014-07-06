function [W,idx_supp,C] = tak_sim_weight_6d_MNI_MTL(clusterPath,q,GRID)
% tak_function
%=========================================================================%
% - Comments
%=========================================================================%
% (07/06/2014)
%%
if ~exist('GRID','var') || isempty(GRID)
    GRID = 'Grid326'; % <- I mostly focus on this...
end

%=========================================================================%
% anomalous edge indices
%=========================================================================%
load(clusterPath,'wsupp', 'idx_anomEdgeClusters', 'idx_anomNodeClusters')

%=========================================================================%
% graph info of the node parcellation
%=========================================================================%
graphPath=[get_rootdir, '/data_local/graphinfo/graph_info_',GRID,'.mat'];
load(graphPath,'C')
%% create weight vector
p = size(C,2);

idx_supp = false(p,1);
idx_supp(wsupp) = true;
% return
%=========================================================================%
% assign ground truth weight on the support
%=========================================================================%
W = zeros(p,q);
%-------------------------------------------------------------------------%
for j=1:q
    % assign weights "cluster-by-cluster"
    for i=1:size(idx_anomEdgeClusters,1)
        idx_tmp = idx_anomEdgeClusters{i,1};

        % make weight piecewise constant + perturbation of noise
        W(idx_tmp,j) = tak_sample_signed_unif([4,10],1) ...
                      + 0*randn([length(idx_tmp),1]);
    end
end