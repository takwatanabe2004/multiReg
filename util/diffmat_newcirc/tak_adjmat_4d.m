function [A,AX,AY,AZ,AT]=tak_adjmat_4d(ARRAYSIZE,flagcirc)
% [A,AX,AY,AZ,AT]=tak_adjmat_4d(ARRAYSIZE,flagcirc)
%-------------------------------------------------------------------------%
% Create adjacency matrix corresonding to the 4-d finite difference matrix
%-------------------------------------------------------------------------%
% INPUT
%   ARRAYSIZE = [X Y Z T] -> dimension size of the array
%   flagcirc: '0' -> make non-circulant difference matrix [default]
%             '1' -> make circulant difference matrix 
%                    (creates wrap-around terms, so use with caution)
%-------------------------------------------------------------------------%
% OUTPUT
%  A: A = AX+AY+AZ+AT
%   where AX,AY,AZ,AT: adjacency graph in the (X,Y,Z,T) direction
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
T=ARRAYSIZE(4);

%=========================================================================%
% Create 1-D adjacency matrix for each dimension
%=========================================================================%
AX1d=tak_adjmat_1d(X,flagcirc);
AY1d=tak_adjmat_1d(Y,flagcirc);
AZ1d=tak_adjmat_1d(Z,flagcirc);
AT1d=tak_adjmat_1d(T,flagcirc);

%=========================================================================%
% create adjacency graph for each array dimension
%-------------------------------------------------------------------------%
% exploit kronecker structure needed to create the difference operator for 
% each dimension (see my research notes)
%=========================================================================%
IX=speye(X);
IY=speye(Y);
IZ=speye(Z);
IT=speye(T);
Iyx=kron(IY,IX);
Itz=kron(IT,IZ);

% create first order difference operator for each array dimension
AX=kron(Itz,kron(IY,AX1d));
AY=kron(Itz,kron(AY1d,IX));
AZ=kron(kron(IT,AZ1d),Iyx);
AT=kron(kron(AT1d,IZ),Iyx);

%=========================================================================%
% create final adjacency graph
%=========================================================================%
A=AX+AY+AZ+AT;