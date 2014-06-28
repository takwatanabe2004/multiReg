function A=tak_adjmat_1d(n,flagcirc)
% A=tak_adjmat_1d(n,flagcirc)
% (06/28/2014)
%==================================================================================
% - create 1D adjacency matrix corresponding to the 1D differencing matrix
%   option=0 
%    -> D = (n-1 x n) matrix with no circulant structure
%   option=1 
%    -> D = (n x n) matrix with circulant structure...
%    -> the circulant structure creates a "wrap-around term"
%       (for an N-D vector x, the last entry of the product D*x will be x1-xN...
%        which may not be sensible)
%==================================================================================
if nargin==1 
    flagcirc=0;
end

% (n-1 x n) incidence/difference matrix 
A=spdiags(ones(n-1,1),-1,n,n); % lower major diagonal part (later to be symmetrized)

if flagcirc % pad a row at the end that adds the (x1-xN) effect
    A(n,1) = 1;
end
A = A+A'; % symmetrize