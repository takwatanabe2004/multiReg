function [X,Z,DISTR] = ...
    tak_sample_connectome3d_batch(NSAMP, T, NSIZE, rt, rspace)
% [X,Z,DISTR] = ...
%     tak_sample_connectome3d_batch(NSAMP, T, rt, rspace, NSIZE)
%=========================================================================%
% "Batches" together tak_prec_connectome2d.m and tak_sample_connectome2d.m
% in a single step.
%-------------------------------------------------------------------------%
% T = # time points
% rt=temporal correlation
% rspace = [rx,ry,rz] <- spatial correlation
% NSIZE = [nx,ny,nz]
%=========================================================================%
% (06/27/2014)
%% create matrix normal distribution
DISTR = tak_prec_connectome3d(T,NSIZE,  rt, rspace);
[X,Z]=tak_sample_connectome3d(DISTR,NSAMP);