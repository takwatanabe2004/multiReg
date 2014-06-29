function [C,Cx,Cy]=tak_diffmat_2d_nokron(NSIZE,flagcirc)
% default: non-circulant difference matrix
if nargin==1 
    flagcirc=0;
end
% NDIM = length(NSIZE);
p = prod(NSIZE);

%% circulant case
if flagcirc 
    %=================================================================%
    % differencing operator in the X-direction
    %=================================================================%
    idx_row = [];
    idx_col = [];
    for i=1:NSIZE(2)
        offset = (i-1)*NSIZE(1);
        % last term is the "wrap-around" part
        idx_row = [idx_row, offset+2:offset+NSIZE(1),   offset+1  ];
        idx_col = [idx_col, offset+1:offset+NSIZE(1)-1, offset+NSIZE(1) ];
%             keyboard
    end
    Cx = speye(p) - sparse(idx_row, idx_col, 1, p,p);

    %=================================================================%
    % differencing operator in the Y-direction
    %=================================================================%
    % 2nd part in these index-set is the "wrap-around" part
    idx_row = [1+NSIZE(1):p,  1:NSIZE(1)]; 
    idx_col = [1:p-NSIZE(1),  p-NSIZE(1)+1:p  ]; 
    Cy = speye(p) - sparse(idx_row, idx_col, 1, p,p);

%% noncirculant case
else 
    %=================================================================%
    % differencing operator in the X-direction
    %=================================================================%
    idx_row = [];
    idx_col = [];
    for i=1:NSIZE(2)
        offset_row = (i-1)*(NSIZE(1)-1);
        offset_col = (i-1)*(NSIZE(1));
        idx_row = [idx_row, offset_row+1:offset_row+NSIZE(1)-1];
        idx_col = [idx_col, offset_col+1:offset_col+NSIZE(1)-1];
%             keyboard
    end
    nrowx = (NSIZE(1)-1) *  NSIZE(2);    
    Cx =  sparse(idx_row, 1+idx_col,  1, nrowx, p) ...
        - sparse(idx_row,   idx_col,  1, nrowx, p);

    %=================================================================%
    % differencing operator in the Y-direction
    %=================================================================%
    nrowy =  NSIZE(1)    * (NSIZE(2)-1);
    idx = 1:nrowy;
    Cy =  sparse(idx, idx+NSIZE(1), 1, nrowy, p) ...
        - sparse(idx, idx,          1, nrowy, p);
end
C = [Cx; Cy];