%% sim_plot_node_2d_improved - loopless version
% (06/26/2014)
%=========================================================================%
% - Plot (nx by ny) node orientation for simulation.
%=========================================================================%
%%
clear all;
purge

nx = 12;
ny = 15;

figure,imexpl,tak_plot_sim_nodes2d(nx,ny)
%% plot anomalous nodes
% first cluster
idx_nx1 = [2:5];
idx_ny1 = [3:6];
[xx1,yy1]=ndgrid(idx_nx1,idx_ny1);
idx_anom1 = sub2ind([nx,ny],xx1(:),yy1(:));

% 2nd cluster
idx_nx2 = [nx-10:nx-8];
idx_ny2 = [ny-5:ny-3];
[xx2,yy2]=ndgrid(idx_nx2,idx_ny2);
idx_anom2 = sub2ind([nx,ny],xx2(:),yy2(:));

idx_anom = [idx_anom1; idx_anom2];

figure,imexpl,tak_plot_sim_nodes2d(nx,ny,idx_anom)
