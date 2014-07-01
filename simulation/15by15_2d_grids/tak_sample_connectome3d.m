function [X,Z] = tak_sample_connectome3d(DISTR, NSAMP)
% [X,Z] = tak_sample_connectome3d(DISTR, NSAMP)
%% (07/01/2014) Created on top of tak_sample_connectome3d.m
%=========================================================================%
% - Sample spatial time series, where nodes are in 3d space.
%   The spatial and temporal correlation are specified in DISTR
%   (3-D AR model for spatial smoothness, 1-D AR model for temporal smoothness)
% - X = (NSAMP x p), where p = # edges
% - Z = (X x Y x Z X T x NSAMP) shaped times series
%-------------------------------------------------------------------------%
% - DISTR: distribution struct created using tak_prec_connectome3d.m
%=========================================================================%
%%
d = DISTR.X*DISTR.Y*DISTR.Z; % # nodes
p = nchoosek(d,2); % # edges
dim = size(DISTR.A,1);



if nargout==2
    Z = zeros(DISTR.T, DISTR.X, DISTR.Y, DISTR.Z,  NSAMP);
    % Z = zeros(DISTR.X, DISTR.Y, DISTR.T, NSAMP);
end

X = zeros(NSAMP,p);
for isamp=1:NSAMP
    tmp=(DISTR.A\randn(dim, 1)); 

    if nargout==2
        Z(:,:,:,:,isamp) = reshape(tmp,[DISTR.T, DISTR.X, DISTR.Y, DISTR.Z]);
%         Z(:,:,:,isamp) = reshape(tmp,[DISTR.X, DISTR.Y, DISTR.T]);
    end

    rcorr_tmp = corr(reshape(tmp,[DISTR.T,d]));
    X(isamp,:) = tak_dvec(rcorr_tmp);
end


%-------------------------------------------------------------------------%
% have the index-configured as (X x Y x Z X T x NSAMP)
if nargout==2
    Z=permute(Z,[2,3,4,1,5]);
end