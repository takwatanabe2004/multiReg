%% june27_6d_adjmat_try.m
% (06/27/2014)
%=========================================================================%
% - trying to create "adjmat" generator for 4d
% - in the same spirit with the "diffmat" generator
% - sanity check via inc2adj, adj2inc
%=========================================================================%
%%
clear all;
purge

nx1=11;
ny1=21;
nz1=13;
nx2=nx1;
ny2=ny1;
nz2=nz1;
NSIZE=[nx1,ny1,nz1,nx2,ny2,nz2];
n=prod(NSIZE);
n/1e6
% C=tak_diffmat(NSIZE,0);
% L=C'*C;
% A=tak_inc2adj(C);

% figure,imexpb
% subplot(131),tspy(C)
% subplot(132),tspy(L)
% subplot(133),tspy(A)
% return
%% try to create adjacency matrix "A" directly
Ax1=spdiags(ones(nx1-1,1),-1,nx1,nx1);
Ax1 = Ax1+Ax1';
Ay1=spdiags(ones(ny1-1,1),-1,ny1,ny1);
Ay1 = Ay1+Ay1';
Az1=spdiags(ones(nz1-1,1),-1,nz1,nz1);
Az1 = Az1+Az1';
Ax2=spdiags(ones(nx2-1,1),-1,nx2,nx2);
Ax2 = Ax2+Ax2';
Ay2=spdiags(ones(ny2-1,1),-1,ny2,ny2);
Ay2 = Ay2+Ay2';
Az2=spdiags(ones(nz2-1,1),-1,nz2,nz2);
Az2 = Az2+Az2';

IX1=speye(nx1);
IY1=speye(ny1);
IZ1=speye(nz1);
IX2=speye(nx2);
IY2=speye(ny2);
IZ2=speye(nz2);

%-------------------------------------------------------------------------%
IZ1_IY1_IX1=kron(IZ1,kron(IY1,IX1));
IZ2_IY2_IX2=kron(IZ2,kron(IY2,IX2));

IZ1_IY1_Ax1=kron(IZ1,kron(IY1,Ax1));
IZ1_Ay1_IX1=kron(IZ1,kron(Ay1,IX1));
Az1_IY1_IX1=kron(Az1,kron(IY1,IX1));

IZ2_IY2_Ax2=kron(IZ2,kron(IY2,Ax2));
IZ2_Ay2_IX2=kron(IZ2,kron(Ay2,IX2));
Az2_IY2_IX2=kron(Az2,kron(IY2,IX2));
%-------------------------------------------------------------------------%
Ax1_blk = kron(IZ2_IY2_IX2,IZ1_IY1_Ax1);
Ay1_blk = kron(IZ2_IY2_IX2,IZ1_Ay1_IX1);
Az1_blk = kron(IZ2_IY2_IX2,Az1_IY1_IX1);
Ax2_blk = kron(IZ2_IY2_Ax2,IZ1_IY1_IX1);
Ay2_blk=kron(IZ2_Ay2_IX2,IZ1_IY1_IX1);
Az2_blk=kron(Az2_IY2_IX2,IZ1_IY1_IX1);
%-------------------------------------------------------------------------%
A2 = Ax1_blk+Ay1_blk+Az1_blk+Ax2_blk+Ay2_blk+Az2_blk;
Whos A A2 C
% figure,imexpb
% subplot(131),tspy(Ax_blk)
% subplot(132),tspy(Ay_blk)
% subplot(133),tspy(A2)
if isequal(A,A2)
    disp('********** success! **********')
else
    error('i suck...')
end
