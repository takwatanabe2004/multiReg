function [W,output]=tak_EN_OCP1_regr_ADMM(X,Y,lam,gam,eta,options,F,Wtrue)
% [W,output]=tak_admm_EN_regr_MTL(X,Y,lam,gam,options,Wtrue)
%=========================================================================%
% - ADMM elastic net regression - MTL:
%    1/2||Y-Xw||^2 + lambda \sum_j||W_j||_2 + gamma/2||W||^2
%                                           + lam \sum_j ||Fw_j||_1
%=========================================================================%
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
%=========================================================================%
% (07/09/2014)
%% sort out 'options'
[n,p]=size(X);
q=size(Y,2);

ef=size(F,1);
if(~exist('options','var')||isempty(options)),     
    rho = 1;
    
    maxiter = 500;
    tol = 5e-4;
    progress = inf;
    silence = false;
    funcval = false;

    if p > n
        K=tak_admm_inv_lemma(X,1/(rho));
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
            K=tak_admm_inv_lemma(X,1/(rho));
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
V2 = zeros(p,q);
V3 = zeros(ef,p);

% dual variables
U1 = zeros(p,q);
U2 = zeros(p,q);
U3 = zeros(ef,p);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
XtY=(X'*Y);
if n >= p
    XtX = X'*X;
end
FtF = F'*F;

% the right-multiplying inverse term for the W-update
M = (gam+2*rho)*speye(q) + rho*FtF; 
% Minv = inv(M);

%=========================================================================%
% keep track of function value
%=========================================================================%
vec = @(W)W(:);
fval = @(W) 1/2 * norm( vec(Y-X*W) )^2 + gam/2*norm(W(:))^2 ...
                                       + lam*norm(W(:),1) ...
                                       + eta*norm(vec(F*W'),1);
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
W_old=W;
% disp('go')

% keep track of function value
if funcval, fvalues=zeros(maxiter,1); end;
if exist('Wtrue','var'), wdist=zeros(maxiter,1); end;
for k=1:maxiter
    if funcval,  fvalues(k)=fval(W); end;   
    if exist('Wtrue','var'), wdist(k)=norm(W(:)-Wtrue(:)); end;
    
    if mod(k,progress)==0 && k~=1
        str='--- %3d out of %d ... Tol=%2.2e (tinner=%4.3fsec, ttotal=%4.3fsec) ---\n';
        fprintf(str,k,maxiter,rel_change,toc(time.inner),toc(time.total))
        time.inner = tic;
    end
    
    %======================================================================
    % update first variable block: (W)
    %======================================================================
%     W = rho*(V1-U1 + V2-U2 + (V3-U3)'*F)*Minv;    
    TMP = (V1-U1 + V2-U2 + (V3-U3)'*F)';
    W = rho*(M\TMP)';
    
    %======================================================================
    % update second variable block: (V1,V2,V3)
    %======================================================================
    % update V1
    Q=(XtY + rho*(W+U1));
    if p > n        
        V1 = Q/(rho) - 1/(rho)^2*(K*(X*Q));
    else
        V1 = (XtX + rho*speye(p))\Q;
    end
    
    % update V2
    V2 = soft(W+U2,lam/rho);
    
    % update V3
    V3 = soft(F*W'+U3,eta/rho);
    
    %======================================================================
    % dual updates
    %======================================================================
    U1=U1+(W-V1);
    U2=U2+(W-V2);
    U3=U3+(F*W'-V3);

    %======================================================================
    % Check termination criteria
    %======================================================================
    %%% relative change in primal variable norm %%%
    rel_change=norm(W(:)-W_old(:))/norm(W_old(:));
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

time.total=toc(time.total);
%% organize output
% primal variables
W=V2;
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
if funcval,  
    fvalues(k+1)=fval(W); % <- final function value
    fvalues=fvalues(1:k+1);
    output.fval=fvalues;
end;

% (optional) distance to wtrue
if exist('Wtrue','var')
    wdist(k+1)=norm(W(:)-Wtrue(:)); 
    wdist=wdist(1:k+1);
    output.wdist=wdist;
end;
