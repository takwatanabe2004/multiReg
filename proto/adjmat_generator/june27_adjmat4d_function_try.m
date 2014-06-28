%% june27_adjmat4d_function_try.m
% (06/27/2014)
%=========================================================================%
% - trying to create "adjmat" generator for 4d
% - in the same spirit with the "diffmat" generator
% - sanity check via inc2adj, adj2inc
% - here use function matlab "function" i created
%=========================================================================%
%%
clear all;
purge

nx=11;
ny=21;
nz=25;
nt=33;
NSIZE=[nx,ny,nz,nt];
n=prod(NSIZE);
% n/1e6
C=tak_diffmat(NSIZE,0);
L=C'*C;
A=tak_inc2adj(C);

Whos C L A
% figure,imexpb
% subplot(131),tspy(C)
% subplot(132),tspy(L)
% subplot(133),tspy(A)
% return
%% try to create adjacency matrix "A" directly
A2 = tak_adjmat(NSIZE,0);
% 
% figure,imexpb
% subplot(131),tspy(Ax_blk)
% subplot(132),tspy(Ay_blk)
% subplot(133),tspy(A2)
if isequal(A,A2)
    disp('********** success! **********')
else
    error('i suck...')
end
%% repeat above for circulant case
Ccirc=tak_diffmat(NSIZE,1);
Lcirc=Ccirc'*Ccirc;
Acirc=tak_inc2adj(Ccirc);
% figure,imexpb
% subplot(131),tspy(Ccirc)
% subplot(132),tspy(Lcirc)
% subplot(133),tspy(Acirc)

Acirc2 = tak_adjmat(NSIZE,1);
if isequal(Acirc,Acirc2)
    disp('********** success! **********')
else
    error('i suck...')
end