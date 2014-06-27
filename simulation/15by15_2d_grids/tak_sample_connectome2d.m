function [X,Z] = tak_sample_connectome2d(DISTR, NSAMP)
% [X,Z] = tak_sample_connectome2d(DISTR)
%=========================================================================%
% - Sample spatial time series, where nodes are in 2d space.
%   The spatial and temporal correlation are specified in DISTR
%   (2-D AR model for spatial smoothness, 1-D AR model for temporal smoothness)
% - X = (p X NSAMP), where p = # edges
% - Z = (X x Y x T x NSAMP) shaped times series
%-------------------------------------------------------------------------%
% - DISTR: distribution struct created using tak_prec_connectome2d.m
%=========================================================================%
% (06/27/2014)
%%
d = DISTR.X*DISTR.Y; % # nodes
p = nchoosek(d,2); % # edges
dim = size(DISTR.A,1);

X = zeros(p,NSAMP);

if nargout==2
    Z = zeros(DISTR.T, DISTR.X, DISTR.Y,  NSAMP);
    % Z = zeros(DISTR.X, DISTR.Y, DISTR.T, NSAMP);
end

for isamp=1:NSAMP
    tmp=(DISTR.A\randn(dim, 1)); 

    if nargout==2
        Z(:,:,:,isamp) = reshape(tmp,[DISTR.T, DISTR.X, DISTR.Y]);
%         Z(:,:,:,isamp) = reshape(tmp,[DISTR.X, DISTR.Y, DISTR.T]);
    end

    rcorr_tmp = corr(reshape(tmp,[DISTR.T,d]));
    X(:,isamp) = tak_dvec(rcorr_tmp);
end


%-------------------------------------------------------------------------%
% have the index-configured as (X x Y x T x NSAMP)
if nargout==2
    Z=permute(Z,[2,3,1,4]);
end