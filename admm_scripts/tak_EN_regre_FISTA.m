function [w,output]=tak_EN_regre_FISTA(X,y,lam,gam,options,wtrue)
% [w,output]=tak_EN_regre_FISTA(X,y,lam,gam,options,wtrue)
% (06/20/2014)
%=========================================================================%
% - ISTA Elastic-Net regression:
%    1/2||y-Xw||^2 + lam * ||w||_1 + gam/2 * ||w||^2
%=========================================================================%
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
%% sort out 'options'
p=size(X,2);

%=========================================================================%
% ISTA paramter and termination criteria
%=========================================================================%
if(~exist('options','var')||isempty(options)),     
    maxiter = 500;
    tol = 5e-4;
    progress = inf;
    silence = false;
    funcval = false;
else
    % step size (needs knowledge of spectral norm of hessian)
    tau = options.tau;

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
w  = zeros(p,1); 
v  = zeros(p,1); 

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
Xty=(X'*y);
Xt=X';

%=========================================================================%
% gradient function handle
%=========================================================================%
GRAD = @(w) Xt*(X*w) - Xty + gam * w;

%=========================================================================%
% keep track of function value (optional, as it could slow down algorithm)
%=========================================================================%
if funcval
    fval = @(w) 1/2 * norm(y-X*w)^2 + lam*norm(w,1) + gam/2*norm(w)^2;
end
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
w_old=w;
t=1;
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
    
    %======================================================================
    % FISTA step
    %=====================================================================%
    w = soft(v - tau*GRAD(v), lam*tau);
    t_old = t;
    t = (1+sqrt(1+4*t^2))/2;
    v = w + ((t_old-1)/t) * (w - w_old);
%     keyboard

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
