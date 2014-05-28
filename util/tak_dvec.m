function dvec = tak_dvec(S)
% - convert the (dxd) symmetric matrix S into d(d-1)/2 vector, extracting 
%   the lower triangular part matrix S.
%%
d = size(S,1);
dvec = S(tril(true(d),-1));
