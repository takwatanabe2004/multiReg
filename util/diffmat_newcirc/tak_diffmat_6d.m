function [C,Cdim]=tak_diffmat_6d(ARRAYSIZE,flagcirc)
% [C,Cdim]=tak_diffmat_6d(ARRAYSIZE,flagcirc)
% comment incomplete...
%----------------------------------------------------------------------------------
% Create 6-d first order difference matrix
%----------------------------------------------------------------------------------
% INPUT
%   ARRAYSIZE = [X1 Y1 Z1 X2 Y2 Z2] -> dimension size of the array
%   flagcirc: '0' -> make non-circulant difference matrix [default]
%             '1' -> make circulant difference matrix 
%                    (creates wrap-around terms, so use with caution)
%----------------------------------------------------------------------------------
% OUTPUT
%  C: C = [CX1; CY1; CZ1; CX2; CY2; CZ2]
%   where CX1,CY1,CZ1,CX2,CY2,CZ2: difference operator on each of the 6 dimensions   
%----------------------------------------------------------------------------------
% Note: if flagcirc==1: 
%   C: (6N x N) matrix, where N = (X1*Y1*Z1*X2*Y2*Z2) 
%   CX1,CY1,CZ1,CX2,CY2,CZ2: (N x N) matrix
%----------------------------------------------------------------------------------
% (05/28/2014)
%----------------------------------------------------------------------------------
%%
% default: non-circulant difference matrix
if nargin==1 
    flagcirc=0;
end

X1=ARRAYSIZE(1);
Y1=ARRAYSIZE(2);
Z1=ARRAYSIZE(3);
X2=ARRAYSIZE(4);
Y2=ARRAYSIZE(5);
Z2=ARRAYSIZE(6);

%==================================================================================
% Create 1-D difference matrix for each dimension
%==================================================================================
DX1=tak_diffmat_1d(X1,flagcirc);
DY1=tak_diffmat_1d(Y1,flagcirc);
DZ1=tak_diffmat_1d(Z1,flagcirc);
DX2=tak_diffmat_1d(X2,flagcirc);
DY2=tak_diffmat_1d(Y2,flagcirc);
DZ2=tak_diffmat_1d(Z2,flagcirc);

%==================================================================================
% create kronecker structure needed to create the difference operator for 
% each dimension (see my research notes)
%==================================================================================
IX1=speye(X1);
IY1=speye(Y1);
IZ1=speye(Z1);
IX2=speye(X2);
IY2=speye(Y2);
IZ2=speye(Z2);

IZ1_IY1_IX1=kron(IZ1,kron(IY1,IX1));
IZ2_IY2_IX2=kron(IZ2,kron(IY2,IX2));

IZ1_IY1_DX1=kron(IZ1,kron(IY1,DX1));
IZ1_DY1_IX1=kron(IZ1,kron(DY1,IX1));
DZ1_IY1_IX1=kron(DZ1,kron(IY1,IX1));

IZ2_IY2_DX2=kron(IZ2,kron(IY2,DX2));
IZ2_DY2_IX2=kron(IZ2,kron(DY2,IX2));
DZ2_IY2_IX2=kron(DZ2,kron(IY2,IX2));

%==================================================================================
% create first order difference operator for each array dimension
%==================================================================================
CX1=kron(IZ2_IY2_IX2,IZ1_IY1_DX1);
CY1=kron(IZ2_IY2_IX2,IZ1_DY1_IX1);
CZ1=kron(IZ2_IY2_IX2,DZ1_IY1_IX1);
CX2=kron(IZ2_IY2_DX2,IZ1_IY1_IX1);
CY2=kron(IZ2_DY2_IX2,IZ1_IY1_IX1);
CZ2=kron(DZ2_IY2_IX2,IZ1_IY1_IX1);

%==================================================================================
% create final difference matrix
%==================================================================================
C=[CX1;CY1;CZ1; CX2; CY2;CZ2];

if nargout ==2
    Cdim.CX1=CX1;
    Cdim.CY1=CY1;
    Cdim.CZ1=CZ1;
    Cdim.CX2=CX2;
    Cdim.CY2=CY2;
    Cdim.CZ2=CZ2;
end 