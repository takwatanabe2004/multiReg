%% sim_plot_node_2d_timeseries_and_edges
% (06/26/2014)
%=========================================================================%
% - Plot time series assigned on the 2-d node orientation, where the
%   nodes have a (nx-ny) orientation
%=========================================================================%
%%
clear all;
purge

nx = 10;
ny = 10;

d = nx*ny; % number of nodes
p = nchoosek(d,2); % number of correlations/edges

%-------------------------------------------------------------------------%
% two clusters of nodes
%-------------------------------------------------------------------------%
% first cluster
idx_nx1 = [8:10]-1;
idx_ny1 = [8:10]-1;

% 2nd cluster
idx_nx2 = [3:5]-1;
idx_ny2 = [3:5]-1;

% figure,imexpl,tak_plot_sim_nodes2d(nx,ny)
% return
%% assign two cluster of nodes
% first cluster
[xx1,yy1]=ndgrid(idx_nx1,idx_ny1);
idx_anom1 = sub2ind([nx,ny],xx1(:),yy1(:))';

% 2nd cluster
[xx2,yy2]=ndgrid(idx_nx2,idx_ny2);
idx_anom2 = sub2ind([nx,ny],xx2(:),yy2(:))';

idx_anom = [idx_anom1, idx_anom2];

figure,imexpl,tak_plot_sim_nodes2d(nx,ny,idx_anom)
% return
%% time series parameters
T  = 200;       % number of time points
rt = 0.0;       % temporal correlation
rx = 0.9;      % spatial correlation (x)
ry = 0.5;      % spatial correlation (y)
nsamp=10;

% % create precision matrix (2d) from the (vectorized) matrix normal distribution
% DISTR = tak_prec_connectome2d(T, [nx,ny],rt,[rx,ry]);
% [X,Z]=tak_sample_connectome2d(DISTR,nsamp);

% sample from precision matrix
[X,Z]=tak_sample_connectome2d_batch(nsamp,T,[nx,ny],rt,[rx,ry]);
Whos Z

%-------------------------------------------------------------------------%
% plot few examples to show spatial & temporal correlation pattern
%-------------------------------------------------------------------------%
tak_plot_neighbortime(Z(:,:,:,1),8,8)
% tak_plot_neighbortime8nn(Z(:,:,:,1),8,8)

imcovvl( tak_dvecinv( X(:,1),1))
% imcovvl( tak_dvecinv( mean(X'),1))
%% display complete bipartite graph of the above two clusters in "connectome space"
mask=false(d,d);
for ii=1:length(idx_anom1)
    idx1 = idx_anom1(ii);
    for jj=1:length(idx_anom2)
        idx2=idx_anom2(jj);
        mask(idx1,idx2)=true;
        mask(idx2,idx1)=true; %<- for symmetry
    end
end
imedger(mask)
%%
[idx,maskVec,maskMat]= tak_nodes2edges_2d(nx,ny,idx_anom1,idx_anom2);
figure,imexpb
subplot(131),tplot(idx)
subplot(132),tplot(maskVec)
subplot(133),imedge(maskMat)

figure,imexpb
subplot(131),tplot(idx)
subplot(132),imedge(tak_dvecinv(maskVec,0))
subplot(133),tplot(tak_dvec(maskMat))