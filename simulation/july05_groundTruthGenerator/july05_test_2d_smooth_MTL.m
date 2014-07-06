%% july05_test_2d_smooth_MTL.m
% (07/05/2014)
%=========================================================================%
% - Test out tak_sim_weight_2d_smooth_MTL.m
%=========================================================================%
%%
clear all;
purge

nx = 256;
ny = 256;
q=20;

[W,idx_supp] = tak_sim_weight_2d_smooth_MTL(nx,ny,q);

figure,imexpl,imagesc(W),impixelinfo,drawnow

figure,imexp
subplot(241),imcov(reshape(W(:,1),[nx,ny])',1)
subplot(242),imcov(reshape(W(:,2),[nx,ny])',1)
subplot(243),imcov(reshape(W(:,3),[nx,ny])',1)
subplot(244),imcov(reshape(W(:,4),[nx,ny])',1)
subplot(245),imcov(reshape(W(:,5),[nx,ny])',1)
subplot(246),imcov(reshape(W(:,6),[nx,ny])',1)
subplot(247),imcov(reshape(W(:,7),[nx,ny])',1)
subplot(248),imcov(reshape(W(:,8),[nx,ny])',1)