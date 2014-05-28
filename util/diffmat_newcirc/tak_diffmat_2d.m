function [C,CX,CY]=tak_diffmat_2d(ARRAYSIZE,flagcirc)
% [C,CX,CY]=tak_diffmat_2d(ARRAYSIZE,flagcirc
%----------------------------------------------------------------------------------
% Create 2-d first order difference matrix
%----------------------------------------------------------------------------------
% INPUT
%   ARRAYSIZE = [X Y Z] -> dimension size of the array
%   flagcirc: '0' -> make non-circulant difference matrix [default]
%             '1' -> make circulant difference matrix 
%                    (creates wrap-around terms, so use with caution)
%----------------------------------------------------------------------------------
% OUTPUT
%  C: C = [CX; CY]
%   where CX,CY: difference operator on each of the 2 dimensions   
%----------------------------------------------------------------------------------
% Note: if flagcirc==1: 
%   C: (2N x N) matrix, where N = (X*Y) 
%   CX,CY: (N x N) matrix
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

%==================================================================================
% Create 1-D difference matrix for each dimension
%==================================================================================
DX=tak_diffmat_1d(X,flagcirc);
DY=tak_diffmat_1d(Y,flagcirc);

%==================================================================================
% create first order difference operator for each array dimension
%==================================================================================
CX=kron(speye(Y),DX);
CY=kron(DY,speye(X));

%==================================================================================
% create final difference matrix
%==================================================================================
C=[CX;CY];