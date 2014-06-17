clear
disp('----------------')

%-------------------------------------------------------------------------%
% p=50e3;
% C=tak_diffmat_1d(p,1);
% LAP  = C'*C + speye(p);
%-------------------------------------------------------------------------%
nx=222;
ny=222;
p=nx*ny;
C = tak_diffmat_2d([nx,ny],1); 
LAP  = C'*C + speye(p);


% [U,D,V]=svds(C,'econ');
% A = sparse(randn(5)+100*eye(5));
% cholmod2(A+A',randn(5,1))

tic
L=lchol(LAP);
toc


tic
LL=chol(LAP')';
toc