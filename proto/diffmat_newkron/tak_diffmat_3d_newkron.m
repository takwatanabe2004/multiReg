function [C,CX,CY,CZ]=tak_diffmat_3d_newkron(ARRAYSIZE,flagcirc)
% [C,CX,CY,CZ]=tak_diffmat_3d_newkron(ARRAYSIZE,flagcirc)
% comment incomplete...
%----------------------------------------------------------------------------------
% Create 3-d first order difference matrix
%----------------------------------------------------------------------------------
% INPUT
%   ARRAYSIZE = [X Y Z] -> dimension size of the array
%   flagcirc: '0' -> make non-circulant difference matrix [default]
%             '1' -> make circulant difference matrix 
%                    (creates wrap-around terms, so use with caution)
%----------------------------------------------------------------------------------
% OUTPUT
%  C: C = [CX; CY; CZ]
%   where CX,CY,CZ: difference operator on each of the 3 dimensions   
%----------------------------------------------------------------------------------
% Note: if flagcirc==1: 
%   C: (3N x N) matrix, where N = (X*Y*Z) 
%   CX,CY,CZ: (N x N) matrix
%----------------------------------------------------------------------------------
% (05/28/2014)
%----------------------------------------------------------------------------------
%%
% default: non-circulant difference matrix
if nargin==1 
    flagcirc=0;
end

X=ARRAYSIZE(1);
Y=ARRAYSIZE(2);
Z=ARRAYSIZE(3);

%==================================================================================
% Create 1-D difference matrix for each dimension
%==================================================================================
DX=tak_diffmat_1d(X,flagcirc);
DY=tak_diffmat_1d(Y,flagcirc);
DZ=tak_diffmat_1d(Z,flagcirc);

%==================================================================================
% create kronecker structure needed to create the difference operator for 
% each dimension (see my research notes)
%==================================================================================
Ix=speye(X);
Iy=speye(Y);
Iz=speye(Z);

%==================================================================================
% create first order difference operator for each array dimension
%==================================================================================
% CX=kron(speye(Z*Y),DX);
% CY=kron(Iz,kron(DY,Ix));
% CZ=kron(DZ,speye(Y*X));
% CX = sparse(KronProd({DX,  1,  1}, [1 2 3], [X,Y,Z],1));
% CY = sparse(KronProd({ 1, DY,  1}, [1 2 3], [X,Y,Z],1));
% CZ = sparse(KronProd({ 1,  1, DZ}, [1 2 3], [X,Y,Z],1));
CX = KronProd({DX,  1,  1}, [1 2 3], [X,Y,Z],1);
CY = KronProd({ 1, DY,  1}, [1 2 3], [X,Y,Z],1);
CZ = KronProd({ 1,  1, DZ}, [1 2 3], [X,Y,Z],1);

% CX=CX.sparse;
% CY=CY.sparse;
% CZ=CZ.sparse;
% keyboard
%==================================================================================
% create final difference matrix
%==================================================================================
% C=[CX;CY;CZ];
C=vertcat(CX,CY,CZ);
p = prod(ARRAYSIZE);
keyboard
DX2 = vertcat(DX, sparse(2*p));
DY2 = vertcat(sparse(p),DY,sparse(p));
DZ2 = vertcat(sparse(2*p),CZ);
keyboard
% C = KronProd({DX,  1,  1}, [1 2 3], [X,Y,Z],1);












