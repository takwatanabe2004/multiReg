function [C,Cx,Cy,Cz,Ct,Cu,Cv]=tak_diffmat_6d_nokron(NSIZE,flagcirc)
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
U = NSIZE(5);
V = NSIZE(6);
XY = X*Y;
YZ = Y*Z;
XYZ = X*Y*Z;
YZT = Y*Z*T;
XYZT=X*Y*Z*T;
XYZTU = X*Y*Z*T*U;
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
    for i=1:prod(NSIZE(2:end))
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
    for i=1:prod(NSIZE(3:end))
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
    for i=1:prod(NSIZE(4:end))
        offset = (i-1)*XYZ;
        idx_row = [idx_row, offset+base_row];
        idx_col = [idx_col, offset+base_col];
    end
    Cz = speye(p) - sparse(idx_row, idx_col, 1, p, p);

    %=================================================================%
    % differencing operator in the T-direction
    %=================================================================%
    idx_row = [];
    idx_col = [];
        
    % 2nd part in these index-set is the "wrap-around" part
    base_row = [XYZ+1 : XYZT   ,              1 : XYZ]; 
    base_col = [    1 : XYZ*(T-1),  XYZ*(T-1)+1 : XYZT];
    for i=1:prod(NSIZE(5:end))
        offset = (i-1)*XYZT;
        idx_row = [idx_row, offset+base_row];
        idx_col = [idx_col, offset+base_col];
    end
    Ct = speye(p) - sparse(idx_row, idx_col, 1, p, p);
    
    %=================================================================%
    % differencing operator in the U-direction
    %=================================================================%
    idx_row = [];
    idx_col = [];
        
    % 2nd part in these index-set is the "wrap-around" part
    base_row = [XYZT+1 : XYZTU   ,              1 : XYZT]; 
    base_col = [    1 : XYZT*(U-1),  XYZT*(U-1)+1 : XYZTU];
    for i=1:prod(NSIZE(6:end))
        offset = (i-1)*XYZTU;
        idx_row = [idx_row, offset+base_row];
        idx_col = [idx_col, offset+base_col];
    end
    Cu = speye(p) - sparse(idx_row, idx_col, 1, p, p);
    
    %=================================================================%
    % differencing operator in the V-direction
    %=================================================================%    
    % 2nd part in these index-set is the "wrap-around" part
    idx_row = [XYZTU+1 : p,                 1 : XYZTU]; 
    idx_col = [    1 : XYZTU*(V-1),  XYZTU*(V-1)+1 : p];
    Cv = speye(p) - sparse(idx_row, idx_col, 1, p, p);

%% noncirculant case
else 
    %=================================================================%
    % differencing operator in the X-direction
    %=================================================================%
    nrowx = (X-1)*Y*Z*T*U*V;    
    idx_row = [];
    idx_col = [];
    base_idx = 1:X-1;
    for i=1:prod(NSIZE(2:end))
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
    nrowy =  X*(Y-1)*Z*T*U*V;
    idx_row = [];
    idx_col = [];    
    base_idx = 1:X*(Y-1);
    for i=1:prod(NSIZE(3:end))
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
    nrowz = X*Y*(Z-1)*T*U*V;
    idx_row = [];
    idx_col = [];    
    base_idx = 1:X*Y*(Z-1);
    for i=1:prod(NSIZE(4:end))
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
    nrowt = X*Y*Z*(T-1)*U*V;
    idx_row = [];
    idx_col = [];    
    base_idx = 1:X*Y*Z*(T-1);
    for i=1:prod(NSIZE(5:end))
        offset_row = (i-1)*XYZ*(T-1);
        offset_col = (i-1)*XYZT;
        idx_row = [idx_row, offset_row+base_idx];
        idx_col = [idx_col, offset_col+base_idx];
    end
    Ct =  sparse(idx_row, idx_col+XYZ, 1, nrowt, p) ...
        - sparse(idx_row, idx_col    , 1, nrowt, p);

    %=================================================================%
    % differencing operator in the T-direction
    %=================================================================%
    nrowu = X*Y*Z*T*(U-1)*V;
    idx_row = [];
    idx_col = [];    
    base_idx = 1:X*Y*Z*T*(U-1);
    for i=1:prod(NSIZE(6:end))
        offset_row = (i-1)*XYZT*(U-1);
        offset_col = (i-1)*XYZTU;
        idx_row = [idx_row, offset_row+base_idx];
        idx_col = [idx_col, offset_col+base_idx];
    end
    Cu =  sparse(idx_row, idx_col+XYZT, 1, nrowu, p) ...
        - sparse(idx_row, idx_col     , 1, nrowu, p);

    %=================================================================%
    % differencing operator in the T-direction
    %=================================================================%
    nrowv = XYZTU*(V-1);
    idx = 1:nrowv;
    Cv =  sparse(idx, idx+XYZTU, 1, nrowv, p) ...
        - sparse(idx, idx   ,  1, nrowv, p);
end
%%
% C=[Cx;Cy;Cz];
% vertcat faster
C=vertcat(Cx,Cy,Cz,Ct,Cu,Cv);