function [W,output]=tak_GN_OCP2_regr_FISTA(X,Y,lam,gam,eta,options,C,F,wtrue)
% [W,output]=tak_GN_OCP2_regr_FISTA(X,Y,lam,gam,eta,options,C,F,wtrue)
%=========================================================================%
% - FISTA elastic net + OCP2 regression:
%    1/2||Y-Xw||^2 + lambda \sum_j||W_j||_2 + gamma/2||W||^2
%                                           + eta/2 ||Fw_j||_1
%=========================================================================%
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
%=========================================================================%
% (07/09/2014)
%% sort out 'options'
p=size(X,2);
q=size(Y,2);

%=========================================================================%
% ISTA paramter and termination criteria
%=========================================================================%
if(~exist('options','var')||isempty(options)),        
    maxiter = 500;
    tol = 5e-4;
    progress = inf;
    silence = false;
    funcval = false;
    
    % step size (needs knowledge of spectral norm of hessian)
    tau = 1/(tnormest(X)^2+gam*tnormest(C'*C)+eta*tnormest(F'*F));
else
    % step size (needs knowledge of spectral norm of hessian)
    if isfield(options,'tau')
        tau = options.tau;
    else
        tau = 1/(tnormest(X)^2+gam*tnormest(C'*C)+eta*tnormest(F'*F));
    end

    %=====================================================================%
    % termination criterion
    %=====================================================================%
    maxiter   = options.maxiter;     % <- maximum number of iterations
    tol       = options.tol;         % <- relative change in the primal variable
    progress  = options.progress;    % <- display "progress" (every k iterations)
    silence   = options.silence;     % <- display termination condition
    funcval   = options.funcval;     % <- track function values (may slow alg.)
end
%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
W  = zeros(p,q); 
V  = zeros(p,q); % <- needed for fista step

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
Xty=(X'*Y);
Xt=X';
FtF=F'*F;
CtC=C'*C;

%=========================================================================%
% gradient function handle
%=========================================================================%
GRAD = @(W) Xt*(X*W) - Xty + gam*(CtC*W) + eta*(W*FtF);

%=========================================================================%
% keep track of function value (optional, as it could slow down algorithm)
%=========================================================================%
vec = @(W)W(:);

fval = @(W) 1/2 * norm( vec(Y-X*W) )^2 + gam/2*norm(vec(C*W))^2 ...
                                       + lam*norm(W(:),1) ...
                                       + eta/2*norm(vec(F*W'))^2;
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
W_old=W;
t=1; % <- needed for fista
% disp('go')

% keep track of function value
if funcval, fvalues=zeros(maxiter,1); end;
if exist('wtrue','var'), wdist=zeros(maxiter,1); end;
for k=1:maxiter
    if funcval,  fvalues(k)=fval(W); end;   
    if exist('wtrue','var'), wdist(k)=norm(W(:)-wtrue(:)); end;
    
    if mod(k,progress)==0 && k~=1
        str='--- %3d out of %d ... Tol=%2.2e (tinner=%4.3fsec, ttotal=%4.3fsec) ---\n';
        fprintf(str,k,maxiter,rel_change,toc(time.inner),toc(time.total))
        time.inner = tic;
    end
    
    %======================================================================
    % FISTA step
    %=====================================================================%
    W = soft(V - tau*GRAD(V),lam*tau);
    t_old = t;
    t = (1+sqrt(1+4*t^2))/2;
    V = W + ((t_old-1)/t) * (W - W_old);
%     keyboard

    %======================================================================
    % Check termination criteria
    %======================================================================
    %%% relative change in primal variable norm %%%
    rel_change=norm(W(:)-W_old(:))/norm(W_old(:));
    rel_changevec(k)=rel_change;
    time.rel_change=tic;
    
    flag1=rel_change<tol;
    if flag1 && (k>30) % allow 30 iterations of burn-in period
        if ~silence
            fprintf('*** Primal var. tolerance reached!!! tol=%6.3e (%d iter, %4.3f sec)\n',rel_change,k,toc(time.total))
        end
        break
    elseif k==maxiter
        if ~silence
            fprintf('*** Max number of iterations reached!!! tol=%6.3e (%d iter, %4.3f sec)\n',rel_change,k,toc(time.total))
        end
    end     
    
    % needed to compute relative change in primal variable
    W_old=W;
end

time.total=toc(time.total);
%% organize output
% primal variables
% output.w=v2;
% output.v=v;

% dual variables
% output.u=u;

% time it took for the algorithm to converge
% output.time=time.total;

% number of iteration it took to converge
% output.k=k;

% final relative change in the primal variable
% output.rel_change=rel_change;

% the "track-record" of the relative change in primal variable
output.rel_changevec=rel_changevec(1:k);

% (optional) final function value
if funcval,  
    fvalues(k+1)=fval(W); 
    fvalues=fvalues(1:k+1);
    output.fval=fvalues;
end;

% (optional) distance to wtrue
if exist('wtrue','var')
    wdist(k+1)=norm(W(:)-wtrue(:)); % <- final function value
    wdist=wdist(1:k+1);
    output.wdist=wdist;
end;
