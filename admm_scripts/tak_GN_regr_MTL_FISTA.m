function [W,output]=tak_GN_regr_MTL_FISTA(X,Y,lam,gam,options,C,wtrue)
% [w,output]=tak_GN_regr_MTL_FISTA(X,y,lam,gam,options,C,wtrue)
% (06/22/2014)
%=========================================================================%
% - ISTA GraphNet regression:
%    1/2||y-Xw||^2 + lam * ||w||_1 + gam/2 * ||C*w||^2
%=========================================================================%
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
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
    tau = 1/(tnormest(X)^2+gam*tnormest(C'*C));
    MTL = true;
else
    %=====================================================================%
    % termination criterion
    %=====================================================================%
    maxiter   = options.maxiter;     % <- maximum number of iterations
    tol       = options.tol;         % <- relative change in the primal variable
    progress  = options.progress;    % <- display "progress" (every k iterations)
    silence   = options.silence;     % <- display termination condition
    funcval   = options.funcval;     % <- track function values (may slow alg.)

    % step size (needs knowledge of spectral norm of hessian)
    if isfield(options,'tau')
        tau = options.tau;
    else
        tau = 1/(tnormest(X)^2+gam*tnormest(C'*C));
    end
    
    if isfield(options,'MTL')
        MTL   = options.MTL;          % <- use L1/L2 penalty? (default: true)
    else
        MTL = true;
    end
end
%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
W  = zeros(p,q); 
V  = zeros(p,q); 

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder
vsoft=@(W,tau) tsoftvec(W,tau); % vector-soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
XtY=(X'*Y);
Xt=X';
CtC=C'*C;

%=========================================================================%
% gradient function handle
%=========================================================================%
GRAD = @(w) Xt*(X*w) - XtY + gam * (CtC*w);

%=========================================================================%
% keep track of function value (optional, as it could slow down algorithm)
%=========================================================================%
vec = @(W)W(:);

if MTL
    % last term sums the 2-norm of the rows
    fval = @(W) 1/2 * norm(Y-X*W)^2 + gam/2*norm( vec(C*W) )^2 + lam*sum(sqrt(sum(W.^2,2))); 
else
    fval = @(W) 1/2 * norm(Y-X*W)^2 + gam/2*norm( vec(C*W) )^2 + lam*norm(W(:),1);
end
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
W_old=W;
t=1;
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
    if MTL
        W = vsoft(V - tau*GRAD(V), lam*tau);
    else
        W = soft(V - tau*GRAD(V), lam*tau);
    end    
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
    if flag1 && (k>80) % allow 30 iterations of burn-in period
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
