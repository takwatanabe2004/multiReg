function [C,Cx,Cy,Cz,Ct]=tak_diffmat_4d_nokron(NSIZE,flagcirc)
% default: non-circulant difference matrix
if nargin==1 
    flagcirc=0;
end
% NDIM = length(NSIZE);
p = prod(NSIZE);

X = NSIZE(1); 
Y = NSIZE(2);
Z = NSIZE(3);
T = NSIZE(4);
XY = X*Y;
YZ = Y*Z;
XYZ = X*Y*Z;
YZT = Y*Z*T;
XYZT=X*Y*Z*T;
%% circulant case
if flagcirc 
    %=================================================================%
    % differencing operator in the X-direction
    %=================================================================%
    idx_row = [];
    idx_col = [];
    
    % 2nd part in these index-set is the "wrap-around" part
    base_row = [2:    X,   1];
    base_col = [1:(X-1),   X];
    for i=1:YZT
        offset = (i-1)*X;
        idx_row = [idx_row, offset+base_row];
        idx_col = [idx_col, offset+base_col];
    end
    Cx = speye(p) - sparse(idx_row, idx_col, 1, p,p);
%     Cx=speye(p);

    %=================================================================%
    % differencing operator in the Y-direction
    %=================================================================%
    idx_row = [];
    idx_col = [];
    
    % 2nd part in these index-set is the "wrap-around" part
    base_row =[X+1 :  XY,             1 : X];
    tmp_col = [  1 : X*(Y-1), X*(Y-1)+1: XY];
    for i=1:T*Z
        offset = (i-1)*XY;        
        idx_row = [idx_row, offset+base_row];
        idx_col = [idx_col, offset+tmp_col];
    end
    Cy = speye(p) - sparse(idx_row, idx_col, 1, p,p);
%     Cy=speye(p);
    
    %=================================================================%
    % differencing operator in the Z-direction
    %=================================================================%
    idx_row = [];
    idx_col = [];
        
    % 2nd part in these index-set is the "wrap-around" part
    base_row = [XY+1 : XYZ   ,             1 : XY]; 
    base_col = [   1 : XY*(Z-1),  XY*(Z-1)+1 : XYZ];
    for i=1:T
        offset = (i-1)*XYZ;
        idx_row = [idx_row, offset+base_row];
        idx_col = [idx_col, offset+base_col];
    end
    Cz = speye(p) - sparse(idx_row, idx_col, 1, p, p);

    %=================================================================%
    % differencing operator in the T-direction
    %=================================================================%
    % 2nd part in these index-set is the "wrap-around" part
    idx_row = [XYZ+1 : XYZT,                 1 : XYZ]; 
    idx_col = [    1 : XYZ*(T-1),  XYZ*(T-1)+1 : XYZT];
    Ct = speye(p) - sparse(idx_row, idx_col, 1, p, p);
%% noncirculant case
else 
    %=================================================================%
    % differencing operator in the X-direction
    %=================================================================%
    nrowx = (X-1)*YZ*T;    
    idx_row = [];
    idx_col = [];
    base_idx = 1:X-1;
    for i=1:Y*Z*T
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
    nrowy =  X*(Y-1)*Z*T;
    idx_row = [];
    idx_col = [];    
    base_idx = 1:X*(Y-1);
    for i=1:Z*T
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
    nrowz = XY*(Z-1)*T;
    idx_row = [];
    idx_col = [];    
    base_idx = 1:XY*(Z-1);
    for i=1:T
        offset_row = (i-1)*XY*(Z-1);
        offset_col = (i-1)*XYZ;
        idx_row = [idx_row, offset_row+base_idx];
        idx_col = [idx_col, offset_col+base_idx];
    end
    Cz =  sparse(idx_row, idx_col+XY, 1, nrowz, p) ...
        - sparse(idx_row, idx_col   , 1, nrowz, p);

    %=================================================================%
    % differencing operator in the T-direction
    %=================================================================%
    nrowt = XYZ*(T-1);
    idx = 1:nrowt;
    Ct =  sparse(idx, idx+XYZ, 1, nrowt, p) ...
        - sparse(idx, idx   ,  1, nrowt, p);
end
%%
% C=[Cx;Cy;Cz];
% vertcat faster
C=vertcat(Cx,Cy,Cz,Ct);