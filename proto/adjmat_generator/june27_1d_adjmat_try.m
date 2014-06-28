%% june27_1d_adjmat_try.m
% (06/27/2014)
%=========================================================================%
% - trying to create "adjmat" generator for 1d
% - in the same spirit with the "diffmat" generator
% - sanity check via inc2adj, adj2inc
%=========================================================================%
%%
clear all;
purge

n = 100;
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
%% try to create adjacency matrix "A" directly
% A=full(A);
% A2=diag(ones(n-1,1),-1);
% A2 = A2+A2';
% if isequal(A,A2)
%     disp('********** success! **********')
% else
%     error('i suck...')
% end
%% try to create adjacency matrix "A" directly using spdiags
A2=spdiags(ones(n-1,1),-1,n,n);
A2 = A2+A2';
if isequal(A,A2)
    disp('********** success! **********')
else
    error('i suck...')
end