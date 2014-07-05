%% july04_pcg_matrix_loopless.m
% (07/04/2014)
%=========================================================================%
% - See if I can looplessly solve (C'C+I)^-1*R in my own PCG implementation.
% - this could help the speed of my pcg-based multioutput FL algorithm.
% - RESULTS: SUCCESS!!! =)
%-------------------------------------------------------------------------%
% - note this is different from july04_cg_matrix_loopless.m, as this script
%   uses the ichol-based preconditioner
%=========================================================================%
%%
clear all;
purge
disp('==================================================================')
rand('state',0)

maxit=100;
numRep = 30;
% n=10e3; 
% nx=50; ny=60; n=[nx,ny];
% nx=30; ny=20; nz=40; n=[nx,ny,nz];
nx=10; ny=15; nz=15; nt=5; n=[nx,ny,nz,nt];
p=prod(n);
q = 20;

C = tak_diffmat(n,0);
A = C'*C + speye(p);

% incomplete cholesky preconditioner
L = ichol(A);
% [W(:,kk),~] = pcg(PCG.A, B(:,kk), PCG.tol, PCG.maxiter, PCG.precond,PCG.precond', W(:,kk));

% Whos

R = rand(p,q);
%% solve for Y=(A)^-1*R via looping over columns of R
disp('---- looper: CG ----')
tic
for irep = 1:numRep
    Y = zeros(p,q); 
    for k=1:q
%         [Y(:,k)] = pcg(A,R(:,k),[],maxit);
        [Y(:,k),~] = pcg(A,R(:,k),[],maxit);
    end
end
timeLoopCG = toc
err1=norm(tak_vec(A*Y - R))/norm(R(:));
% norm(R(:))
if err1 > 1e-3,    error('meh...'), end;


disp('---- looper: PCG ----')
tic
for irep = 1:numRep
    Y = zeros(p,q); 
    for k=1:q
%         [Y(:,k)] = pcg(A,R(:,k),[],maxit,L,L');
        [Y(:,k),~] = pcg(A,R(:,k),[],maxit,L,L');
    end
end
timeLoopPCG = toc
err1=norm(tak_vec(A*Y - R))/norm(R(:));
% norm(R(:))
if err1 > 1e-3,    error('meh...'), end;
%% solve via kron approach
disp('---- kron: CG ----')
tic
for irep = 1:numRep
    [Y2,~] = pcg(kron(speye(q),A),R(:),[],maxit);
%     [Y2] = pcg(kron(speye(q),A),R(:),[],maxit);
    Y2=reshape(Y2,[p,q]);
end
timeKronCG = toc
err2 = norm(tak_vec(A*Y2-R))/norm(R(:));
if err2 > 1e-3,    error('meh...'), end;



disp('---- kron: PCG ----')
kronL = kron(speye(q),L);
tic
for irep = 1:numRep
    [Y2,~] = pcg(kron(speye(q),A),R(:),[],maxit,kronL,kronL');
%     [Y2] = pcg(kron(speye(q),A),R(:),[],maxit,kronL,kronL');
    Y2=reshape(Y2,[p,q]);
end
timeKronPCG = toc
err2 = norm(tak_vec(A*Y2-R))/norm(R(:));
if err2 > 1e-3,    error('meh...'), end;
% return
%% solve using my pcg
disp('---- My CG ----')
tic
for irep = 1:numRep
    [Y3,iter]=tak_cg_matrix(A,R);
%     iter
end
timeMY_CG=toc

err3 = norm(tak_vec(A*Y3-R))/norm(R(:));
if err3 > 1e-3,    error('meh...'), end;
%%
disp('---- My PCG ----')
LLt=L*L';
tic
for irep = 1:numRep
%     [Y3,iter]=tak_pcg_matrix(A,R,[]);
    [Y3,iter]=tak_pcg_matrix(A,R,L);
%     iter
%     [Y3,iter]=tak_pcg_matrix(A,R,LLt);
end
timeMY_PCG=toc

err3 = norm(tak_vec(A*Y3-R))/norm(R(:));
if err3 > 1e-3,    error('meh...'), end;