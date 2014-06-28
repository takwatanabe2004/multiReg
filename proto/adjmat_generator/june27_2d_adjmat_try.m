%% june27_2d_adjmat_try.m
% (06/27/2014)
%=========================================================================%
% - trying to create "adjmat" generator for 2d
% - in the same spirit with the "diffmat" generator
% - sanity check via inc2adj, adj2inc
%=========================================================================%
%%
clear all;
purge

nx=202;
ny=400;
NSIZE=[nx,ny];
n=prod(NSIZE);
C=tak_diffmat(NSIZE,0);
L=C'*C;
A=tak_inc2adj(C);
% figure,imexpb
% subplot(131),tspy(C)
% subplot(132),tspy(L)
% subplot(133),tspy(A)
%% try to create adjacency matrix "A" directly
Ax=spdiags(ones(nx-1,1),-1,nx,nx);
Ax = Ax+Ax';
Ay=spdiags(ones(ny-1,1),-1,ny,ny);
Ay = Ay+Ay';

Ax_blk = kron(speye(ny),Ax);
Ay_blk = kron(Ay,speye(nx));
A2 = Ax_blk+Ay_blk;
% figure,imexpb
% subplot(131),tspy(Ax_blk)
% subplot(132),tspy(Ay_blk)
% subplot(133),tspy(A2)
if isequal(A,A2)
    disp('********** success! **********')
else
    error('i suck...')
end