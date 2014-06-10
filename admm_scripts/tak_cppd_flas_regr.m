function output=tak_cppd_flas_regr(X,y,options,C,wtrue)
% output=tak_cppd_flas_regr(X,y,options,C,wtrue)
% (06/08/2014)
%=========================================================================%
% - ADMM fused lasso regression
%    1/2||y-Xw||^2 + lambda||w||_1 + gamma||C*w||_1
%=========================================================================%
% options.lambda
% options.gamma
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
%% sort out 'options'
% penalty parameters
lambda=options.lambda;
gamma=options.gamma;

% CPPD parameter
tau=options.tau;
sigma=options.sigma;

%==========================================================================
% termination criterion
%==========================================================================
maxiter   = options.termin.maxiter;     % <- maximum number of iterations
tol       = options.termin.tol;         % <- relative change in the primal variable
progress  = options.termin.progress;    % <- display "progress" (every k iterations)
silence   = options.termin.silence;     % <- display termination condition

%==========================================================================
% Matrix K for inversion lemma (optionally precomputed...saves time during gridsearch)
%==========================================================================
if isfield(options,'K')
    K=options.K;
else
    K=tak_admm_inv_lemma(X,tau);
end
%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
p = size(X,2);
F = [lambda*speye(p);gamma*C];
s=size(F,1);

% primal variable
w =zeros(p,1); 
wbar = zeros(p,1);

% dual variable
u=zeros(s,1);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
tau_Xty=tau*(X'*y);
Ft = F';

%=========================================================================%
% keep track of function value (consider tracking SPP-function value)
%=========================================================================%
fval = @(w) 1/2*norm(y-X*w)^2 + norm(F*w,1);
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
w_old=w;
Fw=F*w;
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
    
%     keyboard
    %======================================================================
    % update dual variable u
    %======================================================================
    % update u
    tmp = u+sigma*F*wbar;
    u = tmp - soft(tmp/sigma,1/sigma);
    
    %======================================================================
    % update primal variable w
    %======================================================================
    r = w - tau*(Ft*u) + tau_Xty;
    w = r - tau*(K*(X*r));

    %======================================================================
    % update wbar
    %======================================================================
    wbar = 2*w - w_old;
    
    %======================================================================
    % Check termination criteria
    %======================================================================
    %%% relative change in primal variable norm %%%
    rel_change=norm(w-w_old)/norm(w_old);
    rel_changevec(k)=rel_change;
    time.rel_change=tic;
    
    flag1=rel_change<tol;
    if flag1 && (k>10) % allow 10 iterations of burn-in period
        if ~silence
            fprintf('*** Primal var. tolerance reached!!! tol=%6.3e (%d iter, %4.3f sec)\n',rel_change,k,toc(time.total))
        end
        break
    end    
    
    % needed to compute relative change in primal variable
    w_old=w;
end
fvalues(k+1)=fval(w); % <- final function value
fvalues=fvalues(1:k+1);
time.total=toc(time.total);
%% organize output
% primal variables
output.w=w;

% dual variables
output.u=u;

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