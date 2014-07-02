function DISTR = tak_prec_connectome3d(T,NSIZE, rtime,rspace)
%% (07/01/2014) Created on top of tak_prec_connectome2d.m
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
%      NSIZE = [X,Y,Z]: size of the 3D image
%   rspace = [rx,ry,rz]: spatial correlation in the X,Y,Z coordinate direction
%----------------------------------------------------------------------------------
% OUTPUT: 
%   - contains matrix normal distribution: MNV(0,AR1D,AR3D)
%   DISTR: struct containing the following fields
%       - ICOV: Inverse covariance of the (vectorized) matrix normal distribution
%       -    A: A^T*A = ICOV
%            (note A\z gives the desired realization, where z = randn(T*X*Y,1))
%       - AR1D: Struct containing fields pertaining to 1D AR distribution 
%         (row cov...U in wikipedia definition of matrix normal)
%               * T, rtime: same as input
%               * ICOV, A: same as DISTR.ICOV & DISTR.A, but for the 1D AR model
%       - AR3D: Struct containing fields pertaining to 3D AR model distribution
%         (col cov...V in wikipedia definition of matrix normal)
%               * rx, ry, rz, X, Y, Z, N: same as input
%               * ICOV, A: same as DISTR.ICOV & DISTR.A, but for the 3D AR model
%----------------------------------------------------------------------------------
% Example:
%  - nothing yet...(07/01/2014)
%----------------------------------------------------------------------------------
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
DISTR.AR3D.R=rspace; % spatial correlation in (X,Y,Z)-direction
DISTR.AR3D.NSIZE=NSIZE;
DISTR.AR3D.N=prod(NSIZE); % array-size
[DISTR.AR3D.ICOV,DISTR.AR3D.A] =tak_prec_ar3d(DISTR.AR3D.NSIZE, DISTR.AR3D.R);

%==================================================================================
% overall inverse covariance of the matrix normal distribution
%==================================================================================
% DISTR.ICOV=kron(DISTR.AR3D.ICOV,DISTR.AR1D.ICOV); % <- i never use this....kron here is computationaly expensive so remove
DISTR.A=kron(DISTR.AR3D.A,DISTR.AR1D.A);

%-------------------------------------------------------------------------%
% New: (06/27/2014)
%-------------------------------------------------------------------------%
DISTR.T = prod(T); % number of time points
DISTR.X = NSIZE(1);
DISTR.Y = NSIZE(2);
DISTR.Z = NSIZE(3);