function [B,BX1,BY1,BX2,BY2]=tak_circmask_4d(ARRAYSIZE)
% [B,BX,BY]=tak_circmask_4d(ARRAYSIZE)
%----------------------------------------------------------------------------------
% Create diagonal binary masking matrix B, which masks the "wrap-around" artifacts 
% from using the circulant difference matrix created from tak_diffmat_4d_ver2(ARRAYSIZE,1)
%----------------------------------------------------------------------------------
% B: (4N x 4N) diagonal binary masking matrix, where N=prod(ARRAYSIZE)
%----------------------------------------------------------------------------------
% (05/28/2014)
%----------------------------------------------------------------------------------
%%
X1=ARRAYSIZE(1);
Y1=ARRAYSIZE(2);
X2=ARRAYSIZE(3);
Y2=ARRAYSIZE(4);
N=prod(ARRAYSIZE);

%==================================================================================
% Create 1-D masking matrix for each dimension
%==================================================================================
maskX1=tak_circmask_1d(X1);
maskY1=tak_circmask_1d(Y1);
maskX2=tak_circmask_1d(X2);
maskY2=tak_circmask_1d(Y2);

IX1=speye(X1);
IY1=speye(Y1);
IX2=speye(X2);
IY2=speye(Y2);

IY1_X1=kron(IY1,IX1);
IY2_X2=kron(IY2,IX2);

BX1=kron(IY2_X2,kron(IY1,maskX1));
BY1=kron(IY2_X2,kron(maskY1,IX1));
BX2=kron(kron(IY2,maskX2),IY1_X1);
BY2=kron(kron(maskY2,IX2), IY1_X1 );

%==================================================================================
% Create final masking matrix 
% - NOTE: blkdiag slow in my experience when dealing with large sparse matrix
%==================================================================================
B=[        BX1, sparse(N,N), sparse(N,N), sparse(N,N); ...
   sparse(N,N),         BY1, sparse(N,N), sparse(N,N); ...
   sparse(N,N), sparse(N,N),         BX2, sparse(N,N); ...
   sparse(N,N), sparse(N,N), sparse(N,N),        BY2];
% B2=blkdiag(BX1,BY1,BX2,BY2);
% if isequal(B,B2), disp('Good!'),end;