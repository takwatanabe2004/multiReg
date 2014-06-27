%% sim_plot_node_2d
% (06/26/2014)
%=========================================================================%
% - Plot (nx by ny) node orientation for simulation.
%=========================================================================%
%%
clear all;
purge

nx = 12;
ny = 15;

% figure,tak_plot_sim_nodes2d(nx,ny)
figure,imexpl,tak_plot_sim_nodes2d(nx,ny)
%% plot anomalous nodes

% idx_nx = [11, 5, 12];
% idx_ny = [3, 8, 10];

%=========================================================================%
% first cluster
%=========================================================================%
idx_nx = [2:5];
idx_ny = [3:6];
cnt=1;
for ix=1:length(idx_nx)
    for iy=1:length(idx_ny)
        idx_nx2(cnt)=idx_nx(ix);
        idx_ny2(cnt)=idx_ny(iy);
        cnt=cnt+1;
    end
end
idx_anom = sub2ind([nx,ny],idx_nx2,idx_ny2);

%=========================================================================%
% 2nd cluster
%=========================================================================%
idx_nx = [nx-10:nx-8];
idx_ny = [ny-5:ny-3];
cnt=1;
for ix=1:length(idx_nx)
    for iy=1:length(idx_ny)
        idx_nx2(cnt)=idx_nx(ix);
        idx_ny2(cnt)=idx_ny(iy);
        cnt=cnt+1;
    end
end

idx_anom = [idx_anom, sub2ind([nx,ny],idx_nx2,idx_ny2)];
% idx=randsample(nx*ny,5);

figure,imexpl,tak_plot_sim_nodes2d(nx,ny,idx_anom)
