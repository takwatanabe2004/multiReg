%% july05_test_1d_smooth.m
% (07/05/2014)
%=========================================================================%
% - Test out tak_sim_weight_1d_smooth.m
%=========================================================================%
%%
clear all;
purge

p = 800;

[w,idx_supp] = tak_sim_weight_1d_smooth(p);
tplott(w)