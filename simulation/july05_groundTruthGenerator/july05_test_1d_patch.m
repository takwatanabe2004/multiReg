%% july05_test_1d_patch.m
% (07/05/2014)
%=========================================================================%
% - Test out tak_sim_weight_1d_patch.m
%=========================================================================%
%%
clear all;
purge

p = 1500;

[w,idx_supp] = tak_sim_weight_1d_patch(p);
tplott(w)