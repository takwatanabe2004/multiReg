function output=tak_admm_elasticnet(X,y,options,K)
% output=tak_admm_elasticnet(X,y,options,K)
% (02/22/2014) - changed name from tak_admm_stl_elasticnet.m
% (12/16/2013) - use 1/n scale loss & primal variable for convergence
%% sort out 'options'
%==========================================================================
% loss function
%==========================================================================
switch options.loss
    case 'hinge1'
        loss=tak_handle_hinge1; % hinge-loss
    case 'hinge2'
        loss=tak_handle_hinge2; % squared-hinge
    case 'hubhinge'
        loss=tak_handle_hinge_huber(options.loss_huber_param); % huberized-hinge
    otherwise
        error('Specified loss is not supported!!!')
end

% penalty parameters
lambda=options.lambda;
gamma=options.gamma;

% augmented lagrangian parameters
rho=options.rho;

%==========================================================================
% termination criterion
%==========================================================================
maxiter   = options.termin.maxiter;     % <- maximum number of iterations
tol       = options.termin.tol;         % <- relative change in the primal variable
progress  = options.termin.progress;    % <- display "progress" (every k iterations)
silence   = options.termin.silence;     % <- display termination condition

%==========================================================================
% Matrix K for inversion lemma
%==========================================================================
if nargin < 4 % compute matrix "K" needed for inversion lemma
    K=tak_admm_inv_lemma(X,rho/(rho+gamma));
end
%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
[n,p]=size(X);

% primal variable
w =zeros(p,1); 
v1=zeros(n,1);
v2=zeros(p,1);

% dual variables
u1=zeros(n,1);
u2=zeros(p,1);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
% YX=diag(y)*X; % bsxfun much faster for this
YX=bsxfun(@times,X,y);
YXt=YX';
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
w_old=w;
% disp('go')
for k=1:maxiter
    if mod(k,progress)==0 && k~=1
        str='--- %3d out of %d ... Tol=%2.2e (tinner=%4.3fsec, ttotal=%4.3fsec) ---\n';
        fprintf(str,k,maxiter,rel_change,toc(time.inner),toc(time.total))
        time.inner = tic;
    end
    
    %======================================================================
    % update first variable block: (w)
    %======================================================================
    % update w
    q=rho*(YXt*(v1-u1))+rho*(v2-u2);
    w=q/(rho+gamma) - (rho/(rho+gamma)^2)*(K*(X*q));
    
    %%% precompute terms used more than once %%%
    YXw=YX*w;

    %======================================================================
    % update second variable block: (v1,v2)
    %======================================================================
    % update v1
    v1 = loss.prox(YXw+u1,1/(n*rho));
        
    % update v2
    v2=soft(w+u2,lambda/rho);
    
    %======================================================================
    % dual updates
    %======================================================================
    u1=u1+(YXw-v1);
    u2=u2+(w-v2);

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
time.total=toc(time.total);
%% organize output
% primal variables
output.w=w;
output.v1=v1;
output.v2=v2;

% dual variables
output.u1=u1;
output.u2=u2;

% time it took for the algorithm to converge
output.time=time.total;

% number of iteration it took to converge
output.k=k;

% final relative change in the primal variable
output.rel_change=rel_change;

% the "track-record" of the relative change in primal variable
output.rel_changevec=rel_changevec(1:k);