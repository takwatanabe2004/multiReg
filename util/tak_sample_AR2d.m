function X = tak_sample_AR2d(nx,ny, rx,ry, nsamp)
% X = tak_sample_AR2d(p, r, nsamp)
%=========================================================================%
% - Generate samples from 1D AR distribution.
%-------------------------------------------------------------------------%
% (nx,ny) = size of feature in (x,y)-directions
% (rx,ry) = spatial smoothness in (x,y)-directions
% nsamp = # samples
%-------------------------------------------------------------------------%
% X = (nsamp x p) design matrix
%=========================================================================%
% (07/07/2014)
%%
[~, A] = tak_prec_ar2d(nx,ny, rx,ry);
z=randn(nx*ny,nsamp);
X=(A\z)';