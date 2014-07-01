function [X,DISTR] = ...
    tak_sample_connectome3d_subsamp_batch(NSAMP, T, NSIZE, rt, rspace,idx_subsamp)
% [X,Z,DISTR] = ...
%     tak_sample_connectome3d_batch(NSAMP, T, rt, rspace, NSIZE)
%=========================================================================%
% - Same as tak_sample_connectome3d_subsamp_batch.m, but sumbample timeseries at
%   edge-coordinates that are within the brain-mask of interest.
%=========================================================================%
%% create matrix normal distribution
DISTR = tak_prec_connectome3d(T,NSIZE,  rt, rspace);
[X]=tak_sample_connectome3d_subsamp(DISTR,NSAMP,idx_subsamp);