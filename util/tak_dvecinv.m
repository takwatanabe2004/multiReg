function S = tak_dvecinv(dvec,option)
% Convert p=d*(d-1)/2 vector into (dxd) symmetric matrix, with 1's filled in the diagonals
% Operation is isomorphic with tak_dvec
%%
% if option==1, set diagonal entries to 1's (default)
% if option==0, set diagonal entries to 0's
%%
% 5/24/2013
% 06/16/2013 -> made it lower triangular
% 06/22/2013 -> added option
%%
if ~exist('option','var')||isempty(option), option=1; end

% solve for p in the equation "nchoosek(d,2) = d*(d-1)/2 = P"
p = length(dvec);
d = (1+sqrt(1+8*p))/2; % d = max(roots([1 -1 -2*p]));

S = tril(ones(d),-1);
S(S==1) = dvec;

% symmetrize
S = S + S';

if option %  fill-in values of ones for the diagonal entries
    S = S+eye(d);
end