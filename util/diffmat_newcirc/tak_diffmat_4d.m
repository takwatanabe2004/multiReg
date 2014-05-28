function [C,CX1,CY1,CX2,CY2]=tak_diffmat_4d(ARRAYSIZE,flagcirc)
% [C,Cdim]=tak_diffmat_4d(ARRAYSIZE,flagcirc)
% comment incomplete...
%----------------------------------------------------------------------------------
% Create 4-d first order difference matrix
%----------------------------------------------------------------------------------
% INPUT
%   ARRAYSIZE = [X1 Y1 X2 Y2] -> dimension size of the array
%   flagcirc: '0' -> make non-circulant difference matrix [default]
%             '1' -> make circulant difference matrix 
%                    (creates wrap-around terms, so use with caution)
%----------------------------------------------------------------------------------
% OUTPUT
%  C: C = [CX1; CY1; CX2; CY2]
%   where CX1,CY1,CX2,CY2: difference operator on each of the 4 dimensions   
%----------------------------------------------------------------------------------
% Note: if flagcirc==1: 
%   C: (4N x N) matrix, where N = (X1*Y1*X2*Y2) 
%   CX1,CY1,CX2,CY2: (N x N) matrix
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
X2=ARRAYSIZE(3);
Y2=ARRAYSIZE(4);

%==================================================================================
% Create 1-D difference matrix for each dimension
%==================================================================================
DX1=tak_diffmat_1d(X1,flagcirc);
DY1=tak_diffmat_1d(Y1,flagcirc);
DX2=tak_diffmat_1d(X2,flagcirc);
DY2=tak_diffmat_1d(Y2,flagcirc);

%==================================================================================
% create kronecker structure needed to create the difference operator for 
% each dimension (see my research notes)
%==================================================================================
IX1=speye(X1);
IY1=speye(Y1);
IX2=speye(X2);
IY2=speye(Y2);

IY1_X1=kron(IY1,IX1);
IY2_X2=kron(IY2,IX2);

%==================================================================================
% create first order difference operator for each array dimension
%==================================================================================
CX1=kron(IY2_X2,kron(IY1,DX1));
CY1=kron(IY2_X2,kron(DY1,IX1));
CX2=kron(kron(IY2,DX2),IY1_X1);
CY2=kron( kron(DY2,IX2), IY1_X1 );

%==================================================================================
% create final difference matrix
%==================================================================================
C=[CX1;CY1; CX2; CY2];