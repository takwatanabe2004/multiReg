% (05/28/2014)
%-------------------------------------------------------------------------%
% this directory contains the finite differencing matrix generator (and masker),
% where the "wrap-around-artifact" of the circulancy is placed on the first row
% of the difference matrix (in the original version, this "wrap-around" 
% occurs at the last row of the finite difference matrix)
%
% This allows efficient implementation of the isotropic TV penalty.
%-------------------------------------------------------------------------%