function [X] = tak_sample_connectome3d_subsamp(DISTR, NSAMP,idx_subsamp)
% [X,Z] = tak_sample_connectome3d(DISTR, NSAMP,idx_subsamp)
%% (07/01/2014) Created on top of tak_sample_connectome3d.m
%=========================================================================%
% - Same as tak_sample_connectome3d.m, but sumbample timeseries at
%   edge-coordinates that are within the brain-mask of interest.
%=========================================================================%
%%
d = DISTR.X*DISTR.Y*DISTR.Z; % # nodes
% p = nchoosek(d,2); % # edges
p = length(idx_subsamp);
dim = size(DISTR.A,1);



% if nargout==2
%     Z = zeros(DISTR.T, DISTR.X, DISTR.Y, DISTR.Z,  NSAMP);
%     % Z = zeros(DISTR.X, DISTR.Y, DISTR.T, NSAMP);
% end

X = zeros(NSAMP,p);
for isamp=1:NSAMP
    tmp=(DISTR.A\randn(dim, 1)); 
%     keyboard

%     if nargout==2
%         Z(:,:,:,:,isamp) = reshape(tmp,[DISTR.T, DISTR.X, DISTR.Y, DISTR.Z]);
% %         Z(:,:,:,isamp) = reshape(tmp,[DISTR.X, DISTR.Y, DISTR.T]);
%     end

    rcorr_tmp = corr(reshape(tmp,[DISTR.T,d]));
%     KK=reshape(1:d^2,[d,d]);
    tmp2 = tak_vec(rcorr_tmp);
    X(isamp,:) = tmp2(idx_subsamp);
end


%-------------------------------------------------------------------------%
% have the index-configured as (X x Y x Z X T x NSAMP)
% if nargout==2
%     Z=permute(Z,[2,3,4,1,5]);
% end