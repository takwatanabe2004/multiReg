function [A,Ax,Ay]=tak_adjmat_2d_nokron(NSIZE,flagcirc)
%  [A,Ax,Ay]=tak_adjmat_2d_nokron(NSIZE,flagcirc)
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
% default: non-circulant difference matrix
if nargin==1 
    flagcirc=0;
end
p = prod(NSIZE);

X = NSIZE(1);
Y = NSIZE(2);
%% circulant case
if flagcirc 

%% noncirculant case
else 
    %=================================================================%
    % adjacency graph in the X-direction
    %=================================================================%
    idx_x = [];
    base_idx = 1:X-1;
    for i=1:Y
        offset = (i-1)*X;
        idx_x = [idx_x, offset+base_idx];
    end
    if nargout >1
        Ax =  sparse([idx_x, idx_x+1], [idx_x+1, idx_x], 1, p, p);
    end

    %=================================================================%
    % adjacency graph in the Y-direction
    %=================================================================%
    idx_y = 1:X*(Y-1);
    if nargout >2
        Ay =  sparse([idx_y, idx_y+X], [idx_y+X, idx_y], 1, p, p);
    end
    
    % final index set
    idx1 = [idx_x, idx_x+1, idx_y, idx_y+X];
    idx2 = [idx_x+1, idx_x, idx_y+X, idx_y];
end
A = sparse(idx1,idx2,1,p,p);