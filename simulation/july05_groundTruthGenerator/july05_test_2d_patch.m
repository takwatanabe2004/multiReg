%% july05_test_2d_patch.m
% (07/05/2014)
%=========================================================================%
% - Test out tak_sim_weight_2d_patch.m
%=========================================================================%
%%
clear all;
purge

nx = 256;
ny = 256;

[w,W,idx_supp] = tak_sim_weight_2d_patch(nx,ny);
imcovvl(W',1)