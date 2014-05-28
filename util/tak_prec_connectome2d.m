function DISTR = tak_prec_connectome2d(T, rtime, NSIZE, rspace)
% DISTR = tak_prec_connectome2d(T, rtime, NSIZE, rspace)
%----------------------------------------------------------------------------------
% Assign the inverse covariance matrix for a "vectorized" version of the 
% matrix normal distribution, where the
%   * columns of the distribution: 2-D spatial signal distribution
%   *    rows of the distribution: temporal components 
%
% The vectorized version of the random variable from this distribution has a 
% covariance described by a kronecker product 
% (hence the inverse-covariance is also described by the kronecker product)
% 
% http://en.wikipedia.org/wiki/Matrix_normal_distribution
%----------------------------------------------------------------------------------
% INPUT
% Temporal covariance parameter (row covariance)
%       T: number of time points
%   rtime: temporal correlation
% 
% Spatial covariance parameter (column covariance)
%      NSIZE = [X,Y]: size of the 2D image
%   rspace = [rx,ry]: spatial correlation in the X and Y coordinate direction
%----------------------------------------------------------------------------------
% OUTPUT: 
%   - contains matrix normal distribution: MNV(0,AR1D,AR2D)
%   DISTR: struct containing the following fields
%       - ICOV: Inverse covariance of the (vectorized) matrix normal distribution
%       -    A: A^T*A = ICOV
%            (note A\z gives the desired realization, where z = randn(T*X*Y,1))
%       - AR1D: Struct containing fields pertaining to 1D AR distribution 
%         (row cov...U in wikipedia definition of matrix normal)
%               * T, rtime: same as input
%               * ICOV, A: same as DISTR.ICOV & DISTR.A, but for the 1D AR model
%       - AR2D: Struct containing fields pertaining to 2D AR model distribution
%         (col cov...V in wikipedia definition of matrix normal)
%               * rx, ry, X, Y, N: same as input
%               * ICOV, A: same as DISTR.ICOV & DISTR.A, but for the 2D AR model
%----------------------------------------------------------------------------------
% Example:
%  run script sanity_check_prec_ar2d.m
%----------------------------------------------------------------------------------
% (12/22/2013) - MAJOR REVISION!!! SEE COMMIT
%   - I divided out so the field contains AR1D, AR2D, ICOV, and A
%%
%==================================================================================
% temporal covariance (row covariance)
%==================================================================================
DISTR.AR1D.T=T;     % time points
DISTR.AR1D.rtime=rtime; % temporal correlation
[DISTR.AR1D.ICOV,DISTR.AR1D.A]=tak_prec_ar1d(DISTR.AR1D.T,DISTR.AR1D.rtime);

%==================================================================================
% spatial covariance (col covariance)
%==================================================================================
DISTR.AR2D.rx=rspace(1); % spatial correlation in X-direction
DISTR.AR2D.ry=rspace(2); % spatial correlation in Y-direction
DISTR.AR2D.X=NSIZE(1);
DISTR.AR2D.Y=NSIZE(2);
DISTR.AR2D.N=NSIZE(1)*NSIZE(2); % array-size
[DISTR.AR2D.ICOV,DISTR.AR2D.A] =tak_prec_ar2d(DISTR.AR2D.X,DISTR.AR2D.Y,DISTR.AR2D.rx,DISTR.AR2D.ry);

%==================================================================================
% overall inverse covariance of the matrix normal distribution
%==================================================================================
DISTR.ICOV=kron(DISTR.AR2D.ICOV,DISTR.AR1D.ICOV);
DISTR.A=kron(DISTR.AR2D.A,DISTR.AR1D.A);