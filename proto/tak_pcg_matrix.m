function [X,iter] = tak_pcg_matrix(A,B,L,tol,maxit,X)
% [X,iter] = tak_cg_matrix(A,B,tol,maxit,X)
%=========================================================================%
% my preconditioned conjugate graidnet algorithm
% (L = preconditioner)
%-------------------------------------------------------------------------%
% - X = solution to AX=B
% - X (input) = intial guess (default: 0 matrix)
%=========================================================================%
% (07/04/2014)
%%
% p = size(A,1);
% q = size(B,2);
[p,q] = size(B);
Bnorm = norm(B(:));
if ~exist('tol','var') || isempty(tol) 
    tol = 1e-6;
end

if ~exist('maxit','var') || isempty(maxit)
    maxit = 500;
end

if ~exist('X','var') || isempty(X)
    X = zeros(p,q);
end

if ~exist('L','var') || isempty(L) 
    warning('Preconditioner not provided: use conjugate gradient')
    L = speye(p);
end
Lt = L';

%=========================================================================%
% initialization
%=========================================================================%
R = B - A*X;
Z = (Lt\(L\R));
% Z = L\R;
P = Z;

%% begin iteration
rz_inProd_old = R(:)'*Z(:);
for k=1:maxit
    %=====================================================================%
    % compute terms used more than once in an iteration
    %=====================================================================%
    AP = A*P;
    alpha = rz_inProd_old/(P(:)'*AP(:));
    X = X + alpha*P;
    R = R - alpha*AP;
    
    %=====================================================================%
    % check termination condition
    %=====================================================================%
    rnormNEW = R(:)'*R(:)  ;
    if sqrt(rnormNEW)/Bnorm < tol
        break
    end       
    Z = (Lt\(L\R));
    rz_inProd_new = R(:)'*Z(:);
    beta = rz_inProd_new/rz_inProd_old;
    P = Z + beta*P;
    rz_inProd_old = rz_inProd_new;
end
iter=k;
% sqrt(rnormNEW)