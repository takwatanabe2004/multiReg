function [w,output]=tak_FL_regr_LADMM_L1_loss(X,y,lam,gam,options,C,wtrue)
% [w,output]=tak_FL_regr_LADMM_L1_loss(X,y,lam,gam,options,C,tau,wtrue)
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
    tau = 1/(1+svds(X,1)^2+normest(C'*C));
    
    maxiter = 500;
    tol = 4e-3;
    progress = 50;
    silence = false;
    funcval = false;
else
    % augmented lagrangian parameters
    rho=options.rho;

    % LADMM stepsize parameter
    if isfield(options,'tau')
        tau = options.tau;
    else
%         tau = 1/(1+svds(X,1)^2+normest(C'*C));
        tau = 1/(1+tnormest(X)^2+normest(C'*C));
        
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

end

%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
% primal variable
w  = zeros(p,1); 
v1 = zeros(n,1);
v2 = zeros(p,1);
v3 = zeros(e,1);

% dual variables
u1 = zeros(n,1);
u2 = zeros(p,1);
u3 = zeros(e,1);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder
% prox_quad = @(t,tau) t/(tau+1) % prox of 1/2 ||u-v||^2_2

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
% Xty=(X'*y);
% if n >= p
%     XtX = X'*X;
% end
Ct=C';
CtC=Ct*C;
Ip=speye(p);
CtC_Ip = CtC+Ip;

%=========================================================================%
% keep track of function value (optional, as it could slow down algorithm)
%=========================================================================%
if funcval
    fval = @(w) norm(y-X*w,1) + lam*norm(w,1) + gam*norm(C*w,1);
%     fval = @(w) 1/2 * norm(y-X*w)^2 + lam*norm(w,1) + gam*norm(C*w,1);
end
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
w_old=w;
Xw = zeros(n,1); % <- used multiple times, so precompute
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
    prox_grad = (X'*(X*w_old + u1-v1)) + CtC_Ip*w_old +u2-v2+Ct*(u3-v3);
%     prox_grad = (X'*(Xw + u1-v1)) + CtC_Ip*w_old +u2-v2+Ct*(u3-v3);
%     prox_grad = X'*(X*w) + (CtC + Ip)*w + X'*(u1-v1) + u2-v2 + Ct*(u3-v3);
%     w = w - tau*prox_grad; 
    w = w_old - tau*prox_grad;

    %======================================================================
    % Update 2nd variable block (v1,v2,v3)
    %======================================================================
    % update v1 (translated soft-threshold
%     Xw = X*w; % <- used multiple times, so precompute
    v1 = y + soft(X*w + u1 - y, 1/rho); % <- prox of shifted L1 loss
%     v1 = (y/rho + Xw+u1)/(1/rho+1); % <- prox of (shifted) L2 loss
    
    % update v2
    v2 = soft(w+u2,lam/rho);
    
    % update v3
    v3 = soft(C*w+u3,gam/rho);
    
    %======================================================================
    % dual updates
    %======================================================================
    u1=u1+(X*w-v1);
    u2=u2+(w-v2);
    u3=u3+(C*w-v3);
%     if mod(k,50)==0,keyboard,end;

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
keyboard
%%
tplott(X*w-v1)
tplott(w-v2)
tplott(C*w-v3)
%%

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
