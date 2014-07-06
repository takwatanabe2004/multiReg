%% july05_test_2d_smooth.m
% (07/05/2014)
%=========================================================================%
% - Test out tak_sim_groundTruthWeight_2d_smooth.m
%=========================================================================%
%%
clear all;
purge

nx = 256;
ny = 256;

[w,W,idx_supp] = tak_sim_weight_2d_smooth(nx,ny);
imcovvl(W',1)