function [P, A] = tak_prec_ar3d( imdim, smoothness)
% [P, A] = tak_prec_ar3d(imdim,smoothness)
%|------------------------------------------------------------------------------|%
%| Description:
%|      Create the precision matrix corresponding to the causal 3d AR(1) model.
%|------------------------------------------------------------------------------|%
%| Input: 
%|      imdim = [nx,ny,nz]
%|      smoothness = [rx,ry,rz]
%|      rx = smoothness parameter in the x-direction
%|      ry = smoothness parameter in the y-direction
%|      ry = smoothness parameter in the z-direction
%|------------------------------------------------------------------------------|%
%| Output: 
%|      P = precision matrix in 'sparse' format.  
%|          Do full(P\speye(N)) to recover covariance matrix.
%|      A = A'*A=P
% http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
%|------------------------------------------------------------------------------|%
%| Modified 2/19/2012
%|      I changed the function input from [nx,ny,nz],[rx,ry,rz] to [imdim], smoothness
%|      for concision/cleaness of code.
%|      
%|      Old: [P, A] = tak_prec_ar3d(nx,ny,nz,rx,ry,rz)
%|      New: [P, A] = tak_prec_ar3d(imdim,smoothness)
%|------------------------------------------------------------------------------|%
%| Created 2/13/2012
%|------------------------------------------------------------------------------|%
nx = imdim(1);
ny = imdim(2);
nz = imdim(3);

rx = smoothness(1);
ry = smoothness(2);
rz = smoothness(3);
if nargout == 1
    Px = tak_prec_ar1d(nx,rx);
    Py = tak_prec_ar1d(ny,ry);
    Pz = tak_prec_ar1d(nz,rz);
    P = kron(Pz, kron(Py,Px));
else
    [Px, Ax] = tak_prec_ar1d(nx,rx);
    [Py, Ay] = tak_prec_ar1d(ny,ry);
    [Pz, Az] = tak_prec_ar1d(nz,rz);
    P = kron(Pz, kron(Py,Px));
    A = kron(Az, kron(Ay,Ax));
end
