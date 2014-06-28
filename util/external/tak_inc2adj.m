function adjmat=tak_inc2adj(incmat)
% adjmat=tak_inc2adj(incmat)
%=========================================================================%
% - my lazy function since inc2adj doesn't return symmetric matrix
%=========================================================================%
% (06/27/2014)
%%
adjmat=inc2adj(incmat);
adjmat=adjmat+adjmat'; % symmetrize