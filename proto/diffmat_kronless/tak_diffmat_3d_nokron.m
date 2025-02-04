function [C,Cx,Cy,Cz]=tak_diffmat_3d_nokron(NSIZE,flagcirc)
% default: non-circulant difference matrix
if nargin==1 
    flagcirc=0;
end
% NDIM = length(NSIZE);
p = prod(NSIZE);

X = NSIZE(1); 
Y = NSIZE(2);
Z = NSIZE(3);
XY = X*Y;
YZ = Y*Z;
%% circulant case
if flagcirc 
    %=================================================================%
    % differencing operator in the X-direction
    %=================================================================%
    idx_row = [];
    idx_col = [];
    
    % 2nd part in these index-set is the "wrap-around" part
    tmp_row = [2:    X,   1];
    tmp_col = [1:(X-1),   X];
    for i=1:YZ
        offset = (i-1)*X;
        idx_row = [idx_row, offset+tmp_row];
        idx_col = [idx_col, offset+tmp_col];
    end
    Cx = speye(p) - sparse(idx_row, idx_col, 1, p,p);
%     Cx=speye(p);

    %=================================================================%
    % differencing operator in the Y-direction
    %=================================================================%
    idx_row = [];
    idx_col = [];
    
    % 2nd part in these index-set is the "wrap-around" part
    tmp_row = [1+X :     XY,            1:X];
    tmp_col = [1   :(X-1)*Y,   (X-1)*Y+1:XY];
    for z=1:Z
        offset = (z-1)*XY;        
        idx_row = [idx_row, offset+tmp_row];
        idx_col = [idx_col, offset+tmp_col];
    end
    Cy = speye(p) - sparse(idx_row, idx_col, 1, p,p);
%     Cy=speye(p);
    
    %=================================================================%
    % differencing operator in the Z-direction
    %=================================================================%
    % 2nd part in these index-set is the "wrap-around" part
    idx_row = [1+XY:p,  1:XY]; 
    idx_col = [1:p-XY,  p-XY+1:p];
    Cz = speye(p) - sparse(idx_row, idx_col, 1, p, p);
%% noncirculant case
else 
    %=================================================================%
    % differencing operator in the X-direction
    %=================================================================%
    nrowx = (X-1)*YZ;    
    idx_row = [];
    idx_col = [];
    base_idx = 1:X-1;
    for i=1:YZ
        offset_row = (i-1)*(X-1);
        offset_col = (i-1)*X;
        idx_row = [idx_row, offset_row+base_idx];
        idx_col = [idx_col, offset_col+base_idx];
    end
    Cx =  sparse(idx_row, 1+idx_col,  1, nrowx, p) ...
        - sparse(idx_row,   idx_col,  1, nrowx, p);

    %=================================================================%
    % differencing operator in the Y-direction
    %=================================================================%
    nrowy =  X*(Y-1)*Z;
    idx_row = [];
    idx_col = [];    
    base_idx = 1:X*(Y-1);
    for i=1:Z
        offset_row = (i-1)*X*(Y-1);
        offset_col = (i-1)*XY;
        idx_row = [idx_row, offset_row+base_idx];
        idx_col = [idx_col, offset_col+base_idx];
    end
    Cy =  sparse(idx_row, idx_col+X, 1, nrowy, p) ...
        - sparse(idx_row, idx_col,   1, nrowy, p);

    %=================================================================%
    % differencing operator in the Z-direction
    %=================================================================%
    nrowz = XY*(Z-1);
    idx = 1:nrowz;
    Cz =  sparse(idx, idx+XY, 1, nrowz, p) ...
        - sparse(idx, idx   , 1, nrowz, p);
end
%%
% C=[Cx;Cy;Cz];
% vertcat faster
C=vertcat(Cx,Cy,Cz);