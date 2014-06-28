%% june27_adjmat3d_function_try.m
% (06/27/2014)
%=========================================================================%
% - trying to create "adjmat" generator for 3d
% - in the same spirit with the "diffmat" generator
% - sanity check via inc2adj, adj2inc
% - here use function matlab "function" i created
%=========================================================================%
%%
clear all;
purge

nx=222;
ny=40;
nz=50;
NSIZE=[nx,ny,nz];
n=prod(NSIZE);
C=tak_diffmat(NSIZE,0);
L=C'*C;
A=tak_inc2adj(C);

Whos C L A
% figure,imexpb
% subplot(131),tspy(C)
% subplot(132),tspy(L)
% subplot(133),tspy(A)
%% try to create adjacency matrix "A" directly
A2=tak_adjmat(NSIZE,0);
% figure,imexpb
% subplot(131),tspy(Ax_blk)
% subplot(132),tspy(Ay_blk)
% subplot(133),tspy(A2)
if isequal(A,A2)
    disp('********** success! **********')
else
    error('i suck...')
end
%% what about "circulant" case?
Ccirc=tak_diffmat(NSIZE,1);
Lcirc=Ccirc'*Ccirc;
Acirc=tak_inc2adj(Ccirc);
% figure,imexpb
% subplot(131),tspy(Ccirc)
% subplot(132),tspy(Lcirc)
% subplot(133),tspy(Acirc)

Acirc2=tak_adjmat(NSIZE,1);
if isequal(Acirc,Acirc2)
    disp('********** success! **********')
else
    error('i suck...')
end