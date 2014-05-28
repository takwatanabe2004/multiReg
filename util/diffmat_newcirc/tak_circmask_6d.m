function [B,Bdim]=tak_circmask_6d(ARRAYSIZE)
% [B,Bdim]=tak_circmask_6d(ARRAYSIZE)
%----------------------------------------------------------------------------------
% Create diagonal binary masking matrix B, which masks the "wrap-around" artifacts 
% from using the circulant difference matrix created from tak_diffmat_6d_ver2(ARRAYSIZE,1)
%----------------------------------------------------------------------------------
% B: (6N x 6N) diagonal binary masking matrix, where N=prod(ARRAYSIZE)
%----------------------------------------------------------------------------------
% (05/28/2014)
%----------------------------------------------------------------------------------
%%
X1=ARRAYSIZE(1);
Y1=ARRAYSIZE(2);
Z1=ARRAYSIZE(3);
X2=ARRAYSIZE(4);
Y2=ARRAYSIZE(5);
Z2=ARRAYSIZE(6);
N=prod(ARRAYSIZE);

%==================================================================================
% Create 1-D masking matrix for each dimension
%==================================================================================
maskX1=tak_circmask_1d(X1);
maskY1=tak_circmask_1d(Y1);
maskZ1=tak_circmask_1d(Z1);
maskX2=tak_circmask_1d(X2);
maskY2=tak_circmask_1d(Y2);
maskZ2=tak_circmask_1d(Z2);

IX1=speye(X1);
IY1=speye(Y1);
IZ1=speye(Z1);
IX2=speye(X2);
IY2=speye(Y2);
IZ2=speye(Z2);

IZ1_IY1_IX1=kron(IZ1,kron(IY1,IX1));
IZ2_IY2_IX2=kron(IZ2,kron(IY2,IX2));

IZ1_IY1_maskX1=kron(IZ1,kron(IY1,maskX1));
IZ1_maskY1_IX1=kron(IZ1,kron(maskY1,IX1));
maskZ1_IY1_IX1=kron(maskZ1,kron(IY1,IX1));

IZ2_IY2_maskX2=kron(IZ2,kron(IY2,maskX2));
IZ2_maskY2_IX2=kron(IZ2,kron(maskY2,IX2));
maskZ2_IY2_IX2=kron(maskZ2,kron(IY2,IX2));

BX1=kron(IZ2_IY2_IX2,IZ1_IY1_maskX1);
BY1=kron(IZ2_IY2_IX2,IZ1_maskY1_IX1);
BZ1=kron(IZ2_IY2_IX2,maskZ1_IY1_IX1);
BX2=kron(IZ2_IY2_maskX2,IZ1_IY1_IX1);
BY2=kron(IZ2_maskY2_IX2,IZ1_IY1_IX1);
BZ2=kron(maskZ2_IY2_IX2,IZ1_IY1_IX1);

%==================================================================================
% Create final masking matrix 
% - NOTE: blkdiag slow in my experience when dealing with large sparse matrix
%==================================================================================
B=[        BX1, sparse(N,N), sparse(N,N), sparse(N,N), sparse(N,N), sparse(N,N); ...
   sparse(N,N),         BY1, sparse(N,N), sparse(N,N), sparse(N,N), sparse(N,N); ...
   sparse(N,N), sparse(N,N),         BZ1, sparse(N,N), sparse(N,N), sparse(N,N); ...
   sparse(N,N), sparse(N,N), sparse(N,N),         BX2, sparse(N,N), sparse(N,N); ...
   sparse(N,N), sparse(N,N), sparse(N,N), sparse(N,N),         BY2, sparse(N,N); ... 
   sparse(N,N), sparse(N,N), sparse(N,N), sparse(N,N), sparse(N,N),        BZ2];

% B2=blkdiag(BX1,BY1,BZ1,BX2,BY2,BZ2);
% if isequal(B,B2), disp('Good!'),end;

if nargout ==2
    Bdim.BX1=BX1;
    Bdim.BY1=BY1;
    Bdim.BZ1=BZ1;
    Bdim.BX2=BX2;
    Bdim.BY2=BY2;
    Bdim.BZ2=BZ2;
end 