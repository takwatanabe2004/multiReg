function [A,AX,AY]=tak_adjmat_2d(ARRAYSIZE,flagcirc)
% [A,AX,AY]=tak_adjmat_2d(ARRAYSIZE,flagcirc)
%-------------------------------------------------------------------------%
% Create adjacency matrix corresonding to the 2-d finite difference matrix
%-------------------------------------------------------------------------%
% INPUT
%   ARRAYSIZE = [X Y] -> dimension size of the array
%   flagcirc: '0' -> make non-circulant difference matrix [default]
%             '1' -> make circulant difference matrix 
%                    (creates wrap-around terms, so use with caution)
%-------------------------------------------------------------------------%
% OUTPUT
%  A: A = AX + AY
%   where AX,AY: adjacency graph in the (X,Y) direction
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

%=========================================================================%
% Create 1-D adjacency matrix for each dimension
%=========================================================================%
AX1d=tak_adjmat_1d(X,flagcirc);
AY1d=tak_adjmat_1d(Y,flagcirc);

%=========================================================================%
% create adjacency graph for each array dimension
%=========================================================================%
AX=kron(speye(Y),AX1d);
AY=kron(AY1d,speye(X));

%=========================================================================%


%=========================================================================%
A=AX+AY;