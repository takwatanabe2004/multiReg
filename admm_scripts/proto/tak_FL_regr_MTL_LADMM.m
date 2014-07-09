function [W,output]=tak_FL_regr_MTL_LADMM(X,Y,lam,gam,options,C,wtrue)
% [w,output]=tak_FL_regr_LADMM_ver2(X,y,lam,gam,options,C,tau,wtrue)
%=========================================================================%
% - same as tak_FL_regr_LADMM.m, but different splitting (see onenote)
%  (specifically, don't split Xw into Xv1)
%=========================================================================%
% - Linaerized-ADMM fused lasso net regression:
%    1/2||y-Xw||^2 + lam * ||w||_1 + gamma2 * ||C*w||_1
%=========================================================================%
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
% options.tau = LADMM stepsize parameter
%-------------------------------------------------------------------------%
% (07/07/2014)
%% sort out 'options'
[n,p]=size(X);
q=size(Y,2);
e=size(C,1);

%=========================================================================%
% AL paramter and termination criteria
%=========================================================================%
if(~exist('options','var')||isempty(options)),     
    % augmented lagrangian parameters
    rho = 1;

    % LADMM stepsize parameter
    tau = 1/(1+tnormest(C'*C));
%     tau = 1/(1+normest(C'*C)^2);
    
    maxiter = 500;
    tol = 4e-3;
    progress = 50;
    silence = false;
    funcval = false;
    MTL = true;

    if p > n
        K=tak_admm_inv_lemma(X,tau/rho);
    end
else
    % augmented lagrangian parameters
    rho=options.rho;

    % LADMM stepsize parameter
    if isfield(options,'tau')
        tau = options.tau;
    else
        tau = 1/(1+tnormest(C'*C));
%         tau = 1/(1+normest(C'*C)^2);
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
    if isfield(options,'MTL')
        MTL   = options.MTL;          % <- use L1/L2 penalty? (default: true)
    else
        MTL = true;
    end

    %=====================================================================%
    % Matrix K for inversion lemma 
    % (optionally precomputed...saves time during gridsearch)
    % (only use inversion lemma when p > n...else solve matrix inverse directly)
    %=====================================================================%
    if p > n
        if isfield(options,'K')
            K=options.K;
        else
            K=tak_admm_inv_lemma(X,tau/rho);
        end
    end
end

%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
% primal variable
W  = zeros(p,q); 
V1 = zeros(p,q);
V2 = zeros(e,q);

% dual variables
U1 = zeros(p,q);
U2 = zeros(e,q);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder
vsoft=@(W,tau) tsoftvec(W,tau); % vector-soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
XtY=(X'*Y);
if n >= p
    XtX = X'*X;
end
Ct=C';
CtC=Ct*C;
Ip=speye(p);
CtC_Ip = CtC+Ip;

%=========================================================================%
% keep track of function value (optional, as it could slow down algorithm)
%=========================================================================%
vec = @(W)W(:);
if MTL
    % last term sums the 2-norm of the rows
    fval = @(W) 1/2 * norm( vec(Y-X*W) )^2 + lam*sum(sqrt(sum(W.^2,2))) + ...
                                           + gam*norm( vec(C*W), 1);
else
    fval = @(W) 1/2 * norm( vec(Y-X*W) )^2 + lam*norm(W(:),1) + ...
                                           + gam*norm( vec(C*W), 1);
end
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
W_old=W;
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
    
    %=====================================================================%
    % update w (linearized ADM step)
    %=====================================================================%
    % L-ADMM prox-gradient term
    prox_grad = CtC_Ip*W_old + (U1-V1)+(Ct*(U2-V2));
    
    % apply inversion lemma
    Q = (tau/rho)*XtY + (W_old - tau*prox_grad);
    if p > n
        W = Q - (tau/rho)*(K*(X*Q));
    else
        W = ( (tau/rho)*XtX + Ip)\Q;
    end

    %======================================================================
    % Update 2nd variable block (v1,v2)
    %======================================================================
    % update v1
    if MTL
        V1 = vsoft(W+U1,lam/rho);
    else
        V1 = soft(W+U1,lam/rho);
    end  
    
    % update v2
    V2 = soft(C*W+U2,gam/rho);
    
    %======================================================================
    % dual updates
    %======================================================================
    U1=U1+(W-V1);
    U2=U2+(C*W-V2);

    %======================================================================
    % Check termination criteria
    %======================================================================
    %%% relative change in primal variable norm %%%
    rel_change=norm(W(:)-W_old(:))/norm(W_old(:));
    rel_changevec(k)=rel_change;
    time.rel_change=tic;
    
    flag1=rel_change<tol;
    if flag1 && (k>100) % allow 30 iterations of burn-in period
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
% keyboard
time.total=toc(time.total);
%% organize output
% primal variables
W=V1;
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
