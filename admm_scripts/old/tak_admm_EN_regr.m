function [w,output]=tak_admm_EN_regr(X,y,lam,gam,options,wtrue)
% [w,output]=tak_admm_EN_regr(X,y,options,wtrue)
% (06/16/2014)
%=========================================================================%
% - ADMM elastic net regression:
%    1/2||y-Xw||^2 + lambda||w||_1 + gamma/2||w||^2
%=========================================================================%
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
%% sort out 'options'
warning('tak_EN_regr_ADMM.m is the state of the art (06/22/2014)')
[n,p]=size(X);
if(~exist('options','var')||isempty(options)),     
    rho = 1;
    
    maxiter = 500;
    tol = 5e-4;
    progress = inf;
    silence = false;

    if p > n
        K=tak_admm_inv_lemma(X,1/(rho+gam));
    end
else
    % augmented lagrangian parameters
    rho=options.rho;

    %=====================================================================%
    % termination criterion
    %=====================================================================%
    maxiter   = options.maxiter;     % <- maximum number of iterations
    tol       = options.tol;         % <- relative change in the primal variable
    progress  = options.progress;    % <- display "progress" (every k iterations)
    silence   = options.silence;     % <- display termination condition

    %=====================================================================%
    % Matrix K for inversion lemma 
    % (optionally precomputed...saves time during gridsearch)
    % (only use inversion lemma when p > n...else solve matrix inverse directly)
    %=====================================================================%
    if p > n
        if isfield(options,'K')
            K=options.K;
        else
            K=tak_admm_inv_lemma(X,1/(rho+gam));
        end
    end
end





%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
% primal variable
w =zeros(p,1); 
v=zeros(p,1);

% dual variables
u=zeros(p,1);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
Xty=(X'*y);
if n >= p
    XtX = X'*X;
end

%=========================================================================%
% keep track of function value
%=========================================================================%
fval = @(w) 1/2 * norm(y-X*w)^2 + gam/2*norm(w)^2 + lam*norm(w,1);
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
w_old=w;
% disp('go')

% keep track of function value
fvalues=zeros(maxiter,1);
if exist('wtrue','var'), wdist=zeros(maxiter,1); end;
for k=1:maxiter
    fvalues(k)=fval(v);
    if exist('wtrue','var'), wdist(k)=norm(w-wtrue); end;
    
    if mod(k,progress)==0 && k~=1
        str='--- %3d out of %d ... Tol=%2.2e (tinner=%4.3fsec, ttotal=%4.3fsec) ---\n';
        fprintf(str,k,maxiter,rel_change,toc(time.inner),toc(time.total))
        time.inner = tic;
    end
    
    %======================================================================
    % update first variable block: (w)
    %======================================================================
    % update w
    if p > n
        q=(Xty + rho*(v-u));
        w=q/(rho+gam) - 1/(rho+gam)^2*(K*(X*q));
    else
        w = (XtX + (rho+gam)*speye(p))\(Xty + rho*(v-u));
    end

    %======================================================================
    % update second variable block: (v)
    %======================================================================
    % update v
    v = soft(w+u,lam/rho);
    
    %======================================================================
    % dual updates
    %======================================================================
    u=u+(w-v);

    %======================================================================
    % Check termination criteria
    %======================================================================
    %%% relative change in primal variable norm %%%
    rel_change=norm(w-w_old)/norm(w_old);
    rel_changevec(k)=rel_change;
    time.rel_change=tic;
    
    flag1=rel_change<tol;
    if flag1 && (k>111) % allow 30 iterations of burn-in period
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
    w_old=w;
end
fvalues(k+1)=fval(v); % <- final function value
fvalues=fvalues(1:k+1);
time.total=toc(time.total);
%% organize output
% primal variables
w=v;
% output.w=v;
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

% function value 
output.fval=fvalues;

% (optional) distance to wtrue
if exist('wtrue','var')
    wdist(k+1)=norm(w-wtrue); % <- final function value
    wdist=wdist(1:k+1);
    output.wdist=wdist;
end;
