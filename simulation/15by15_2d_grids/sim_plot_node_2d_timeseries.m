%% sim_plot_node_2d_timeseries
% (06/26/2014)
%=========================================================================%
% - Plot time series assigned on the 2-d node orientation, where the
%   nodes have a (nx-ny) orientation
%=========================================================================%
%%
clear all;
purge

nx = 12;
ny = 15;

d = nx*ny; % number of nodes
p = nchoosek(d,2); % number of correlations/edges

% figure,tak_plot_sim_nodes2d(nx,ny)
% figure,imexpl,tak_plot_sim_nodes2d(nx,ny)
%% time series parameters
T  = 200;       % number of time points
rt = 0.9;       % temporal correlation
rx = 0.9;      % spatial correlation (x)
ry = 0.8;      % spatial correlation (y)

% create precision matrix (2d) from the (vectorized) matrix normal distribution
DISTR = tak_prec_connectome2d(T, [nx,ny],rt, [rx,ry]);

nsamp=10;
[X,Z]=tak_sample_connectome2d(DISTR,nsamp);
Whos Z

tak_plot_neighbortime(Z(:,:,:,1),8,8)
% tak_plot_neighbortime8nn(Z(:,:,:,1),8,8)

imcovvl( tak_dvecinv( X(:,1),1))
imcovvl( tak_dvecinv( mean(X'),1))