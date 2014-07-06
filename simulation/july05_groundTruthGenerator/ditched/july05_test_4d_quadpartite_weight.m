%% july05_test_4d_quadpartite_weight.m
% (07/05/2014)
%=========================================================================%
% - Test out tak_sim_weight_4d_quadpartite.m
%=========================================================================%
%%
clear all;
purge

nx=15;
ny=15;
d=nx*ny;
p=nchoosek(d,2);
[w,idx_supp,supp_info,C] = tak_sim_weight_4d_quadpartite(nx,ny);
% return
%=========================================================================%
% visualziation
%=========================================================================%
figure,imexpl,tak_plot_sim_nodes2d(nx,ny,supp_info.idx_anom)
% imedgel(supp_info.maskMat)
% imedger(tak_dvecinv(supp_info.mask,0))
% 
figure,imexpb
subplot(131),tplot(idx_supp)
subplot(132),tplot(supp_info.mask)
subplot(133),imconnEdge(supp_info.mask)
drawnow

imconnl(w)