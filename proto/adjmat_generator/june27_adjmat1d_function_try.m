%% june27_adjmat1d_function_try.m
% (06/27/2014)
%=========================================================================%
% - trying to create "adjmat" generator for 1d
% - in the same spirit with the "diffmat" generator
% - sanity check via inc2adj, adj2inc
% - here use function matlab "function" i created
%=========================================================================%
%%
clear all;
purge

n = 99;
C=tak_diffmat_1d(n,0);
L=C'*C;
A=tak_inc2adj(C);
% figure,imexpb
% subplot(131),imcov(C)
% subplot(132),imcov(L)
% subplot(133),imcov(A)
% figure,imexpb
% subplot(131),tspy(C)
% subplot(132),tspy(L)
% subplot(133),tspy(A)
%% try to create adjacency matrix "A" directly using spdiags
A2 = tak_adjmat(n);
if isequal(A,A2)
    disp('********** success! **********')
else
    error('i suck...')
end
%% what about "circulant" case?
Ccirc=tak_diffmat_1d(n,1);
Lcirc=Ccirc'*Ccirc;
Acirc=tak_inc2adj(Ccirc);

Acirc2=tak_adjmat(n,1);
if isequal(Acirc,Acirc2)
    disp('********** success! **********')
else
    error('i suck...')
end
% figure,imexpb
% subplot(131),imcov(Ccirc)
% subplot(132),imcov(Lcirc)
% subplot(133),imcov(Acirc)
% tspyl(A)
% tspyl(Acirc)