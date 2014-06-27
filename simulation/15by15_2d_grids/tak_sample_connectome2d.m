function Z = tak_sample_connectome2d(DISTR, NSAMP)
% Z = tak_sample_connectome2d(DISTR)
%=========================================================================%
% - Sample spatial time series, where nodes are in 2d space.
%   The spatial and temporal correlation are specified in DISTR
%   (2-D AR model for spatial smoothness, 1-D AR model for temporal smoothness)
% - Z = (X x Y x T x NSAMP) shaped times series
%-------------------------------------------------------------------------%
% - DISTR: distribution struct created using tak_prec_connectome2d.m
%=========================================================================%
% (06/27/2014)
%%
% d = DISTR.X*DISTR.Y; % # nodes
dim = size(DISTR.A,1);

Z = zeros(DISTR.X, DISTR.Y, DISTR.T, NSAMP);

for isamp=1:NSAMP
    tmp=(DISTR.A\randn(dim, 1)); 
    Z(:,:,:,isamp) = reshape(tmp,[DISTR.X, DISTR.Y, DISTR.T]);
end
