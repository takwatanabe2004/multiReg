function [B,BX,BY]=tak_circmask_2d(ARRAYSIZE)
% [B,BX,BY]=tak_circmask_2d(ARRAYSIZE)
%----------------------------------------------------------------------------------
% Create diagonal binary masking matrix B, which masks the "wrap-around" artifacts 
% from using the circulant difference matrix created from tak_diffmat_2d_ver2(ARRAYSIZE,1)
%----------------------------------------------------------------------------------
% B: (2N x 2N) diagonal binary masking matrix, where N=prod(ARRAYSIZE)
%----------------------------------------------------------------------------------
% (05/28/2014)
%----------------------------------------------------------------------------------
%%
X=ARRAYSIZE(1);
Y=ARRAYSIZE(2);
N=prod(ARRAYSIZE);

%==================================================================================
% Create 1-D masking matrix for each dimension
%==================================================================================
maskX=tak_circmask_1d(X);
maskY=tak_circmask_1d(Y);

BX=kron(speye(Y),maskX);
BY=kron(maskY,speye(X));

%==================================================================================
% Create final masking matrix 
% - NOTE: blkdiag slow in my experience when dealing with large sparse matrix
%==================================================================================
B=[BX, sparse(N,N); sparse(N,N), BY];

% B2=blkdiag(BX,BY);
% if isequal(B,B2), disp('Good!'),end;