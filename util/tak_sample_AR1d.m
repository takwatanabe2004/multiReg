function X = tak_sample_AR1d(p, r, nsamp)
% X = tak_sample_AR1d(p, r, nsamp)
%=========================================================================%
% - Generate samples from 1D AR distribution.
%-------------------------------------------------------------------------%
% p = size of feature
% r = temporal correlation
% nsamp = # samples
%-------------------------------------------------------------------------%
% X = (nsamp x p) design matrix
%=========================================================================%
% (07/07/2014)
%%
[~, A] = tak_prec_ar1d(p, r);
z=randn(p,nsamp);
X=(A\z)';