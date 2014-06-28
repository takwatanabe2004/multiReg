%% june27_4d_adjmat_try.m
% (06/27/2014)
%=========================================================================%
% - trying to create "adjmat" generator for 4d
% - in the same spirit with the "diffmat" generator
% - sanity check via inc2adj, adj2inc
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

% figure,imexpb
% subplot(131),tspy(C)
% subplot(132),tspy(L)
% subplot(133),tspy(A)
% return
%% try to create adjacency matrix "A" directly
Ax=spdiags(ones(nx-1,1),-1,nx,nx);
Ax = Ax+Ax';
Ay=spdiags(ones(ny-1,1),-1,ny,ny);
Ay = Ay+Ay';
Az=spdiags(ones(nz-1,1),-1,nz,nz);
Az = Az+Az';
At=spdiags(ones(nt-1,1),-1,nt,nt);
At = At+At';

Ix=speye(nx);
Iy=speye(ny);
Iz=speye(nz);
It=speye(nt);

Iyx=kron(Iy,Ix);
Itz=kron(It,Iz);

Ax_blk = kron(Itz,kron(Iy,Ax));
Ay_blk = kron(Itz,kron(Ay,Ix));
Az_blk = kron(kron(It,Az),Iyx);
At_blk = kron(kron(At,Iz),Iyx);
A2 = Ax_blk+Ay_blk+Az_blk+At_blk;
Whos
% figure,imexpb
% subplot(131),tspy(Ax_blk)
% subplot(132),tspy(Ay_blk)
% subplot(133),tspy(A2)
if isequal(A,A2)
    disp('********** success! **********')
else
    error('i suck...')
end
Whos A A2 C