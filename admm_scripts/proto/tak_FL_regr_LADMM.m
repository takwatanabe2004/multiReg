function [w,output]=tak_FL_regr_LADMM(X,y,lam,gam,options,C,wtrue)
% [w,output]=tak_FL_regr_LADMM(X,y,lam,gam,options,C,tau,wtrue)
% (06/16/2014)
%=========================================================================%
% - Linaerized-ADMM fused lasso net regression:
%    1/2||y-Xw||^2 + lam * ||w||_1 + gamma2 * ||C*w||_1
%=========================================================================%
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
% options.tau = LADMM stepsize parameter
%-------------------------------------------------------------------------%
% (07/06/2014)- added PCG.precond option
%% sort out 'options'
[n,p]=size(X);
e=size(C,1);

%=========================================================================%
% AL paramter and termination criteria
%=========================================================================%
if(~exist('options','var')||isempty(options)),     
    % augmented lagrangian parameters
    rho = 1;

    % LADMM stepsize parameter
    tau = 1/(2+normest(C'*C));
%     tau = 1/(2+normest(C'*C)^2);
    
    maxiter = 500;
    tol = 4e-3;
    progress = 50;
    silence = false;
    funcval = false;

    if p > n
        K=tak_admm_inv_lemma(X,1/rho);
    end
else
    % augmented lagrangian parameters
    rho=options.rho;

    % LADMM stepsize parameter
    if isfield(options,'tau')
        tau = options.tau;
    else
        tau = 1/(2+normest(C'*C));
%         tau = 1/(2+normest(C'*C)^2);
        % svds(A,1)^2
        % eigs(A'*A,1);
    end

    %=====================================================================%
    % termination criterion
    %=====================================================================%
    maxiter   = options.maxiter;     % <- maximum number of iterations
    tol       = options.tol;         % <- relative change in the primal variable
    progress  = options.progress;    % <- display "progress" (every k iterations)
    silence   = options.silence;     % <- display termination condition
    funcval   = options.funcval;     % <- track function values (may slow alg.)

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
CtC_2Ip = CtC+2*Ip;

%=========================================================================%
% keep track of function value (optional, as it could slow down algorithm)
%=========================================================================%
if funcval
    fval = @(w) 1/2 * norm(y-X*w)^2 + lam*norm(w,1) + gam*norm(C*w,1);
end
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
w_old=w;
% disp('go')

% keep track of function value
if funcval, fvalues=zeros(maxiter,1); end;
if exist('wtrue','var'), wdist=zeros(maxiter,1); end;
for k=1:maxiter
    if funcval,  fvalues(k)=fval(w); end;   
    if exist('wtrue','var'), wdist(k)=norm(w-wtrue); end;
    
    if mod(k,progress)==0 && k~=1
        str='--- %3d out of %d ... Tol=%2.2e (tinner=%4.3fsec, ttotal=%4.3fsec) ---\n';
        fprintf(str,k,maxiter,rel_change,toc(time.inner),toc(time.total))
        time.inner = tic;
    end    
    
    %=====================================================================%
    % update w (linearized ADM step)
    %=====================================================================%
    % L-ADMM prox-gradient term
    prox_grad = CtC_2Ip*w_old + u1-v1+u2-v2+Ct*(u3-v3);
    w = w_old - tau*prox_grad; 

    %======================================================================
    % Update 2nd variable block (v1,v2,v3)
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

time.total=toc(time.total);
%% organize output
% primal variables
w=v2;
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
    fvalues(k+1)=fval(w); 
    fvalues=fvalues(1:k+1);
    output.fval=fvalues;
end;

% (optional) distance to wtrue
if exist('wtrue','var')
    wdist(k+1)=norm(w-wtrue); % <- final function value
    wdist=wdist(1:k+1);
    output.wdist=wdist;
end;
