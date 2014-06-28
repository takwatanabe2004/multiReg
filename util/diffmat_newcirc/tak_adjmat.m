function A = tak_adjmat(ARRAYSIZE,flagcirc)
% A = tak_adjmat(ARRAYSIZE,flagcirc)
% (06/28/2014)
%--------------------------------------------------------------------------
% A wrapper for making adjacency matrix for n-d tensor signal
% - the circulant matrix has the "wrap-around-effect" occuring on the 
%   first row (the previous version had this on the last row).
% (default: non-circulant case)
%--------------------------------------------------------------------------
%%
if nargin==1 
    flagcirc=0; % (default: non-circulant)
end

switch length(ARRAYSIZE)
    case 1
        A=tak_adjmat_1d(ARRAYSIZE,flagcirc);
    case 2
        A=tak_adjmat_2d(ARRAYSIZE,flagcirc);
    case 3
        A=tak_adjmat_3d(ARRAYSIZE,flagcirc);
    case 4
        A=tak_adjmat_4d(ARRAYSIZE,flagcirc);
    case 6
        A=tak_adjmat_6d(ARRAYSIZE,flagcirc);
    otherwise
        error('Unsupported dimension!!!')
end