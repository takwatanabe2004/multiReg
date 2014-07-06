%% july05_test_4d_pentapartite_weight_MTL.m
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
q=10;
[W,idx_supp,supp_info,C] = tak_sim_weight_4d_pentapartite_MTL(nx,ny,q);
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

figure,imexpl,imagesc(W),impixelinfo,drawnow

figure,imexp
subplot(241),imconn(W(:,1),1)
subplot(242),imconn(W(:,2),1)
subplot(243),imconn(W(:,3),1)
subplot(244),imconn(W(:,4),1)
subplot(245),imconn(W(:,5),1)
subplot(246),imconn(W(:,6),1)
subplot(247),imconn(W(:,7),1)
subplot(248),imconn(W(:,8),1)
% imconnl(W,1)