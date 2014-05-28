function [P, A] = tak_prec_ar2d(nx,ny,rx,ry)
% [P, A] = tak_prec_ar2d(nx,ny,rx,ry)
%|-----------------------------------------------------------------------|%
%| Description:
%|      Create the precision matrix corresponding to the causal 2d AR(1) model.
%|-----------------------------------------------------------------------|%
%| Input: 
%|      nx
%|      ny
%|      rx = smoothness parameter in the x-direction
%|      ry = smoothness parameter in the y-direction
%|-----------------------------------------------------------------------|%
%| Output: 
%|      P = precision matrix in 'sparse' format.  
%|          Do full(P\speye(N)) to recover covariance matrix.
%|      A = A'*A=P
%|-----------------------------------------------------------------------|%
%| http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
%| To sample from matrix normal distribution, simply do this
% z=randn(nx*ny,1);
% x=A\z;
%|-----------------------------------------------------------------------|%
%| Created 2/13/2012
%|-----------------------------------------------------------------------|%
if nargout == 1
    Px = tak_prec_ar1d(nx,rx);
    Py = tak_prec_ar1d(ny,ry);
    P = kron(Py,Px);
else
    [Px, Ax] = tak_prec_ar1d(nx,rx);
    [Py, Ay] = tak_prec_ar1d(ny,ry);
    P = kron(Py,Px);
    A = kron(Ay,Ax);
end
