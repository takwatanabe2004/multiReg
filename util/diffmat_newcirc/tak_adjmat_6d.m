function [A,AX1,AY1,AZ1,AX2,AY2,AZ2]=tak_adjmat_6d(ARRAYSIZE,flagcirc)
% [A,AX1,AY1,AZ1,AX2,AY2,AZ2]=tak_adjmat_6d(ARRAYSIZE,flagcirc)
%-------------------------------------------------------------------------%
% Create adjacency matrix corresonding to the 4-d finite difference matrix
%-------------------------------------------------------------------------%
% INPUT
%   ARRAYSIZE = [X1 Y1 Z1 X2 Y2 Z2] -> dimension size of the array
%   flagcirc: '0' -> make non-circulant difference matrix [default]
%             '1' -> make circulant difference matrix 
%                    (creates wrap-around terms, so use with caution)
%-------------------------------------------------------------------------%
% OUTPUT
%  A: A = AX1+AY1+AZ1+AX2+AY2+AZ2
%   where AX,AY,AZ,AT: adjacency graph in the 6 dimensions
%-------------------------------------------------------------------------%
% (06/28/2014)
%-------------------------------------------------------------------------%
%%
% default: non-circulant adjacency matrix
if nargin==1 
    flagcirc=0;
end

X1=ARRAYSIZE(1);
Y1=ARRAYSIZE(2);
Z1=ARRAYSIZE(3);
X2=ARRAYSIZE(4);
Y2=ARRAYSIZE(5);
Z2=ARRAYSIZE(6);

%=========================================================================%
% Create 1-D adjacency matrix for each dimension
%=========================================================================%
AX1_1d=tak_adjmat_1d(X1,flagcirc);
AY1_1d=tak_adjmat_1d(Y1,flagcirc);
AZ1_1d=tak_adjmat_1d(Z1,flagcirc);
AX2_1d=tak_adjmat_1d(X2,flagcirc);
AY2_1d=tak_adjmat_1d(Y2,flagcirc);
AZ2_1d=tak_adjmat_1d(Z2,flagcirc);

%=========================================================================%
% create adjacency graph for each array dimension
%-------------------------------------------------------------------------%
% exploit kronecker structure needed to create the difference operator for 
% each dimension (see my research notes)
%=========================================================================%
IX1=speye(X1);
IY1=speye(Y1);
IZ1=speye(Z1);
IX2=speye(X2);
IY2=speye(Y2);
IZ2=speye(Z2);

%-------------------------------------------------------------------------%
IZ1_IY1_IX1=kron(IZ1,kron(IY1,IX1));
IZ2_IY2_IX2=kron(IZ2,kron(IY2,IX2));

%-------------------------------------------------------------------------%
IZ1_IY1_AX1=kron(IZ1,kron(IY1,AX1_1d));
IZ1_AY1_IX1=kron(IZ1,kron(AY1_1d,IX1));
AZ1_IY1_IX1=kron(AZ1_1d,kron(IY1,IX1));

IZ2_IY2_AX2=kron(IZ2,kron(IY2,AX2_1d));
IZ2_AY2_IX2=kron(IZ2,kron(AY2_1d,IX2));
AZ2_IY2_IX2=kron(AZ2_1d,kron(IY2,IX2));

%-------------------------------------------------------------------------%
AX1 = kron(IZ2_IY2_IX2,IZ1_IY1_AX1);
AY1 = kron(IZ2_IY2_IX2,IZ1_AY1_IX1);
AZ1 = kron(IZ2_IY2_IX2,AZ1_IY1_IX1);
AX2 = kron(IZ2_IY2_AX2,IZ1_IY1_IX1);
AY2 = kron(IZ2_AY2_IX2,IZ1_IY1_IX1);
AZ2 = kron(AZ2_IY2_IX2,IZ1_IY1_IX1);

%=========================================================================%
% create final adjacency graph
%=========================================================================%
A=AX1+AY1+AZ1+AX2+AY2+AZ2;