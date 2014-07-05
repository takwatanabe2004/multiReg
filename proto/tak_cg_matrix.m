function [X,iter,resid] = tak_cg_matrix(A,B,tol,maxit,X)
% [X,iter] = tak_cg_matrix(A,B,tol,maxit,X)
%=========================================================================%
% my conjugate graidnet algorithm
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

%=========================================================================%
% initialization
%=========================================================================%
R = B - A*X;
P = R;

%% begin iteration

for k=1:maxit
    %=====================================================================%
    % compute terms used more than once in an iteration
    %=====================================================================%
    rnorm = R(:)'*R(:); %
    AP = A*P;
    alpha = rnorm/(P(:)'*AP(:));
    X = X + alpha*P;
%     Rold = R;
    R = R - alpha*AP;
    
    %=====================================================================%
    % check termination condition
    %=====================================================================%
    rnormNEW = R(:)'*R(:);    
%     errr=sqrt(rnormNEW)/Bnorm
    if sqrt(rnormNEW)/Bnorm < tol
%         if ~silence
%             fprintf('*** Tolerance reached!!! tol=%6.3e (%d iter)\n',sqrt(rnormNEW),k)
%         end
        break
    elseif k==maxit
%         if ~silence
%             fprintf('*** Max number of iterations reached!!! tol=%6.3e (%d iter)\n',sqrt(rnormNEW),k)
%         end
    end       
    beta = rnormNEW/rnorm;
    P = R + beta*P;
end
iter=k;
resid=sqrt(rnormNEW)/Bnorm;
% sqrt(rnormNEW)/Bnorm