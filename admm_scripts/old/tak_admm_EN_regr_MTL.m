function [W,output]=tak_admm_EN_regr_MTL(X,Y,lam,gam,options,Wtrue)
% [W,output]=tak_admm_EN_regr_MTL(X,Y,lam,gam,options,Wtrue)
% (06/22/2014)
%=========================================================================%
% - ADMM elastic net regression - MTL:
%    1/2||Y-Xw||^2 + lambda \sum_j||W_j||_2 + gamma/2||W||^2
%=========================================================================%
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
%% sort out 'options'
warning('tak_EN_regr_MTL_ADMM.m is the state of the art (06/22/2014)')
[n,p]=size(X);
q=size(Y,2);
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
W = zeros(p,q); 
V = zeros(p,q);

% dual variables
U=zeros(p,q);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder
vsoft=@(w,tau) tak_softvec(w, tau);

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
XtY=(X'*Y);
if n >= p
    XtX = X'*X;
end

%=========================================================================%
% keep track of function value
%=========================================================================%
vec = @(W)W(:);
% last term sums the 2-norm of the rows
fval = @(W) 1/2 * norm( vec(Y-X*W) )^2 + gam/2*norm(W(:))^2 + lam*sum(sqrt(sum(W.^2,2))); 
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
W_old=W;
% disp('go')

% keep track of function value
fvalues=zeros(maxiter,1);
if exist('Wtrue','var'), wdist=zeros(maxiter,1); end;
for k=1:maxiter
    fvalues(k)=fval(V);
    if exist('Wtrue','var'), wdist(k)=norm(V(:)-Wtrue(:)); end;
    
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
        Q=(XtY + rho*(V-U));
        W=Q/(rho+gam) - 1/(rho+gam)^2*(K*(X*Q));
    else
        W = (XtX + (rho+gam)*speye(p))\(XtY + rho*(V-U));
    end

    %======================================================================
    % update second variable block: (v)
    %======================================================================
    % update v (note: vector soft-threshold faster in my c-implementation)    
%     V = soft(W+U,lam/rho);
%     for j=1:p
%         V(j,:)=vsoft(W(j,:)+U(j,:),lam/rho);
%     end
    V = tsoftvec(W+U,lam/rho);
    
    %======================================================================
    % dual updates
    %======================================================================
    U=U+(W-V);

    %======================================================================
    % Check termination criteria
    %======================================================================
    %%% relative change in primal variable norm %%%
    rel_change=norm(W-W_old)/norm(W_old);
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
    W_old=W;
end
fvalues(k+1)=fval(V); % <- final function value
fvalues=fvalues(1:k+1);
time.total=toc(time.total);
%% organize output
% primal variables
W=V;
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
if exist('Wtrue','var')
    wdist(k+1)=norm(V(:)-Wtrue(:)); 
    wdist=wdist(1:k+1);
    output.wdist=wdist;
end;
