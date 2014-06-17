function output=tak_admm_FL_regr_pcg(X,y,lam,gam,options,C,PCG,wtrue)
% output=tak_admm_FL_regr_pcg(X,y,lam,gam,options,C,PCG,wtrue)
% (06/16/2014)
%=========================================================================%
% - ADMM fused lasso net regression:
%    1/2||y-Xw||^2 + lam * ||w||_1 + gamma2 * ||C*w||_1
%=========================================================================%
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
%% sort out 'options'
[n,p]=size(X);
e=size(C,1);

%=========================================================================%
% AL paramter and termination criteria
%=========================================================================%
if(~exist('options','var')||isempty(options)),     
    rho = 1e-2;
    
    maxiter = 500;
    tol = 5e-4;
    progress = inf;
    silence = false;

    if p > n
        K=tak_admm_inv_lemma(X,1/rho);
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
            K=tak_admm_inv_lemma(X,1/rho);
        end
    end
end

%=========================================================================%
% conjugate gradient parameters
%=========================================================================%
if(~exist('PCG','var')||isempty(PCG))
    PCG.tol = 1e-8;
    PCG.maxiter = 500;
end
%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
% primal variable
w  = zeros(p,1); 
v1 = zeros(p,1);
v2 = zeros(p,1);
v3 = zeros(e,1);

% dual variables
u1 = zeros(p,1);
u2 = zeros(p,1);
u3 = zeros(e,1);

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
Ct=C';
CtC=Ct*C;
Ip=speye(p);
PCG.A = CtC+2*Ip;

%=========================================================================%
% keep track of function value
%=========================================================================%
fval = @(w) 1/2 * norm(y-X*w)^2 + gam/2*norm(w)^2 + lam*norm(C*w,1);
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
    fvalues(k)=fval(w);
    if exist('wtrue','var'), wdist(k)=norm(w-wtrue); end;
    
    if mod(k,progress)==0 && k~=1
        str='--- %3d out of %d ... Tol=%2.2e (tinner=%4.3fsec, ttotal=%4.3fsec) ---\n';
        fprintf(str,k,maxiter,rel_change,toc(time.inner),toc(time.total))
        time.inner = tic;
    end
    
    
    % update w (conjugate gradient)
    b = (v1+v2-u1-u2) + Ct*(v3-u3);
    [w,~] = pcg(PCG.A, b, PCG.tol, PCG.maxiter, [],[], w);
%     if mod(k,20)==0, keyboard, end;

    %======================================================================
    % Update primal variables
    %======================================================================
    % update v1 (if p > n, apply inversion lemma)
    q=Xty + rho*(w+u1);
    if p > n
        v1=q/rho - 1/rho^2*(K*(X*q));
    else
        v1 = (XtX + rho*Ip)\q;
    end
    
    % update v2
    v2 = soft(w+u2,lam/rho);
    
    % update v3
    v3 = soft(C*w+u3,gam/rho);

    
    %======================================================================
    % dual updates
    %======================================================================
    u1=u1+(w-v1);
    u2=u2+(w-v2);
    u3=u3+(C*w-v3);

    %======================================================================
    % Check termination criteria
    %======================================================================
    %%% relative change in primal variable norm %%%
    rel_change=norm(w-w_old)/norm(w_old);
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
    w_old=w;
end
fvalues(k+1)=fval(w); % <- final function value
fvalues=fvalues(1:k+1);
time.total=toc(time.total);
%% organize output
% primal variables
output.w=v2;
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
