function [A,AX,AY,AZ]=tak_adjmat_3d(ARRAYSIZE,flagcirc)
% [A,AX,AY,AZ]=tak_adjmat_3d(ARRAYSIZE,flagcirc)
%-------------------------------------------------------------------------%
% Create adjacency matrix corresonding to the 3-d finite difference matrix
%-------------------------------------------------------------------------%
% INPUT
%   ARRAYSIZE = [X Y, Z] -> dimension size of the array
%   flagcirc: '0' -> make non-circulant difference matrix [default]
%             '1' -> make circulant difference matrix 
%                    (creates wrap-around terms, so use with caution)
%-------------------------------------------------------------------------%
% OUTPUT
%  A: A = AX+AY+AZ
%   where AX,AY,AZ: adjacency graph in the (X,Y,Z) direction
%-------------------------------------------------------------------------%
% (06/28/2014)
%-------------------------------------------------------------------------%
%%
% default: non-circulant adjacency matrix
if nargin==1 
    flagcirc=0;
end

X=ARRAYSIZE(1);
Y=ARRAYSIZE(2);
Z=ARRAYSIZE(3);

%=========================================================================%
% Create 1-D adjacency matrix for each dimension
%=========================================================================%
AX1d=tak_adjmat_1d(X,flagcirc);
AY1d=tak_adjmat_1d(Y,flagcirc);
AZ1d=tak_adjmat_1d(Z,flagcirc);

%=========================================================================%
% create adjacency graph for each array dimension
%-------------------------------------------------------------------------%
% exploit kronecker structure needed to create the difference operator for 
% each dimension (see my research notes)
%=========================================================================%
Ix=speye(X);
Iy=speye(Y);
Iz=speye(Z);
Iyx=kron(Iy,Ix);
Izy=kron(Iz,Iy);

% create first order difference operator for each array dimension
AX=kron(Izy,AX1d);
AY=kron(Iz,kron(AY1d,Ix));
AZ=kron(AZ1d,Iyx);

%=========================================================================%
% create final adjacency graph
%=========================================================================%
A=AX+AY+AZ;