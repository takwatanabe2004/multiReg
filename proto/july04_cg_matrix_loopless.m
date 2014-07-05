%% july04_cg_matrix_loopless.m
% (07/04/2014)
%=========================================================================%
% - See if I can looplessly solve (C'C+I)^-1*R in my own CG implementation.
% - this could help the speed of my cg-based multioutput FL algorithm.
%=========================================================================%
%%
clear all;
purge

randn('state',0)

numRep = 50;
% n=10e3; 
% nx=50; ny=60; n=[nx,ny];
nx=30; ny=20; nz=20; n=[nx,ny,nz];
p = prod(n);
q = 50;

C = tak_diffmat(n,0);
A = C'*C + speye(p);
% Whos

R = randn(p,q);

maxit=100;
%% solve for Y=(A)^-1*R via looping over columns of R
tic
for irep = 1:numRep
    Y = zeros(p,q); 
    for k=1:q
        [Y(:,k),~] = pcg(A,R(:,k),[],maxit);
    end
end
timeLoop = toc
err1 = norm(tak_vec(A*Y-R))/norm(R(:))
% norm(R(:))
if err1 > 1e-3,    error('meh...'), end;
%% solve via kron approach
tic
for irep = 1:numRep
    [Y2,~,relres,iter] = pcg(kron(speye(q),A),R(:),[],maxit);
    Y2=reshape(Y2,[p,q]);
%     iter
%     err2 = norm(tak_vec(A*Y2-R))/norm(R(:));
%     [relres-err2]
end
timeKron = toc
err2 = norm(tak_vec(A*Y2-R))/norm(R(:))
if err2 > 1e-3,    error('meh...'), end;
%% solve using my pcg
tic
for irep = 1:numRep
    [Y3,iter]=tak_cg_matrix(A,R);
%     iter
end
timeMYCG=toc

err3 = norm(tak_vec(A*Y3-R))/norm(R(:))
if err3 > 1e-3,    error('meh...'), end;
    
