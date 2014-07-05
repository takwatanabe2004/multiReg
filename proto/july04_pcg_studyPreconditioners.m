%% july04_pcg_studyPreconditioners.m
% (07/04/2014)
%=========================================================================%
% - Study different ichol based preconditioners
%-------------------------------------------------------------------------%
% - Conlcusion: the default ichol seems the best...other preconditioneres
%   may take less iterations, but the preconditioning matrix gets denser,
%   hence the cost per pcg iteration gets larger for them...
%=========================================================================%
%%
clear all;
purge
disp('==================================================================')
rand('state',0)

maxit=100;
numRep = 5;
% n=10e3; 
% nx=50; ny=60; n=[nx,ny];
% nx=30; ny=20; nz=40; n=[nx,ny,nz];
nx=10; ny=15; nz=15; nt=15; n=[nx,ny,nz,nt];
p=prod(n);
q = 20;

B = rand(p,q);

C = tak_diffmat(n,0);
A = C'*C + speye(p);
%% solve using my cg
disp('---- My CG ----')
tic
SUM = 0;
for irep = 1:numRep
    [Y1,iter1]=tak_cg_matrix(A,B);
    SUM =SUM +iter1;
end
SUM =SUM/numRep
timeMY_CG=toc

err1 = norm(tak_vec(A*Y1-B))/norm(B(:));
if err1 > 1e-6,    error('meh...'), end;
%%
% incomplete cholesky preconditioner: default
L = ichol(A);
disp('---- My PCG ----')
tic
SUM = 0;
for irep = 1:numRep
    [Y2,iter2]=tak_pcg_matrix(A,B,L);
    SUM =SUM +iter2;
end
SUM =SUM/numRep
timeMY_PCG=toc

err2 = norm(tak_vec(A*Y2-B))/norm(B(:));
if err2 > 1e-6,    error('meh...'), end;

%%
L5 = ichol(A, struct('type','ict','droptol',1e-3));
disp('---- My PCG ----')
tic
SUM = 0;
for irep = 1:numRep
    [Y3,iter3]=tak_pcg_matrix(A,B,L5);
    SUM =SUM +iter3;
end
SUM =SUM/numRep
timeMY_PCG3=toc

err3 = norm(tak_vec(A*Y3-B))/norm(B(:));
if err3 > 1e-6,    error('meh...'), end;
%%
alpha =.005;
L_diagcomp = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
disp('---- My PCG ----')
tic
SUM = 0;
for irep = 1:numRep
    [Y4,iter4]=tak_pcg_matrix(A,B,L_diagcomp);
    SUM =SUM +iter4;
end
SUM =SUM/numRep
timeMY_PCG4=toc