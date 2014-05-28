function DISTR = tak_prec_connectome1d(T, rtime, N, rspace)
% DISTR = tak_prec_connectome1d(T, rtime, N, rspace)
%----------------------------------------------------------------------------------
% Assign the inverse covariance matrix for a "vectorized" version of the 
% matrix normal distribution, where the
%   * columns of the distribution: 1-D spatial signal distribution
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
%   T: number of time points
%   rtime: temporal correlation
% 
% Spatial covariance parameter (column covariance)
%   N: length of the 1D image
%   rspace: spatial correlation
%----------------------------------------------------------------------------------
% OUTPUT
%   DISTR: struct containing the following fields
%       - T, rtime, rspace, N: same as input
%       - U, Usqrt: temporal (row) inverse-covariance and its matrix square root,
%                   whose entries are known in closed form 
%                   (note this is not the cholesky...matrix-sqrt is not unique)
%       - V, Vsqrt: spatial (col) inverse-covariance and its matrix square root,
%                   whose entries are known in closed form 
%       - ICOV,ICOVsqrt: final inverse covariance matrix of the distribution 
%                        and its matrix square root
%       - ICOVsqrtTrans: transpote of the matrix square root
% 
% Note:
%  - sampling from the inverse covariance requires to solve x=inv(icov_sqrt')*z,
%    where the matlab "backslash" can be used: x=(ICOVsqrtTrans\z)
%  - the sampling scheme from the matrix normal distribution is more detailed in 
%    t_j10_1d_seed_matrixNormal.m
%----------------------------------------------------------------------------------
% Example:
%
%----------------------------------------------------------------------------------
% 06/19/2013
% 06/23/2013 -> added 'ARRAYSIZE' field
%%
%==================================================================================
% temporal covariance (row covariance)
%==================================================================================
DISTR.T=T;     % time points
DISTR.rtime=rtime; % temporal correlation
[DISTR.U,tmp]=tak_prec_ar1d(DISTR.T,DISTR.rtime);
DISTR.Usqrt=tmp';

%==================================================================================
% spatial covariance (col covariance)
%==================================================================================
DISTR.rspace=rspace; % spatial correlation (1-d space)
DISTR.N=N; % array-size
DISTR.ARRAYSIZE=[DISTR.N]; %<-06/23/2013...useful for code compatibility with 2d script
[DISTR.V,tmp] =tak_prec_ar1d(DISTR.N,DISTR.rspace);
DISTR.Vsqrt=tmp';

%==================================================================================
% overall inverse covariance of the matrix normal distribution
%==================================================================================
DISTR.ICOV=kron(DISTR.V,DISTR.U);
DISTR.ICOVsqrt=kron(DISTR.Vsqrt,DISTR.Usqrt);

% sampling from the inverse covariance requires to solve x=inv(icov_sqrt')*z...ie the transpose is needed
DISTR.ICOVsqrtTrans=DISTR.ICOVsqrt'; 