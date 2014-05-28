function C = tak_diffmat(ARRAYSIZE,flagcirc)
% C = tak_diffmat(ARRAYSIZE,flagcirc)
% (05/28/2014)
%--------------------------------------------------------------------------
% A wrapper for making difference matrix for n-d tensor signal
% - the circulant matrix has the "wrap-around-effect" occuring on the 
%   first row (the previous version had this on the last row).
%--------------------------------------------------------------------------
%%
switch length(ARRAYSIZE)
    case 1
        C=tak_diffmat_1d(ARRAYSIZE,flagcirc);
    case 2
        C=tak_diffmat_2d(ARRAYSIZE,flagcirc);
    case 3
        C=tak_diffmat_3d(ARRAYSIZE,flagcirc);
    case 4
        C=tak_diffmat_4d(ARRAYSIZE,flagcirc);
    case 6
        C=tak_diffmat_6d(ARRAYSIZE,flagcirc);
    otherwise
        error('Unsupported dimension!!!')
end