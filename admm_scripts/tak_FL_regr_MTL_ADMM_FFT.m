function [W,output]=tak_FL_regr_MTL_ADMM_FFT(X,Y,lam,gam,options,graphInfo,wtrue)
% [W,output]=tak_FL_regr_MTL_ADMM_FFT(X,Y,lam,gam,options,graphInfo,wtrue)
% (06/29/2014)
%=========================================================================%
% - ADMM fused lasso net regression:
%    1/2||y-Xw||^2 + lam * ||w||_1 + gamma2 * ||C*w||_1
%=========================================================================%
% options.K <- optionally precompute
% wtrue <- optional...measure norm(west-wtrue) over iterations if inputted
%% sort out 'options'
[n,p]=size(X);
q=size(Y,2);
C = graphInfo.C;
A = graphInfo.A;
b = graphInfo.b; 
B = repmat(b,[1,q]);
NSIZE = graphInfo.NSIZE; 

pp=size(A,1);
e=size(C,1);

%=========================================================================%
% AL paramter and termination criteria
%=========================================================================%
if(~exist('options','var')||isempty(options)),     
    rho = 1;
    
    maxiter = 500;
    tol = 5e-4;
    progress = inf;
    silence = false;
    funcval = false;
    MTL = true;

    if p > n
        K=tak_admm_inv_lemma(X,1/(2*rho));
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
            K=tak_admm_inv_lemma(X,1/(2*rho));
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
V3 = zeros(pp,q);

% dual variables
U1 = zeros(p,q);
U2 = zeros(e,q);
U3 = zeros(pp,q);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder
vsoft=@(W,tau) tsoftvec(W,tau); % vector-soft-thresholder
% bsoft=@(t,tau) soft(t,tau).*b + (~b).*t; % for v3 update

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
Xty=(X'*Y);
if n >= p
    XtX = X'*X;
end
Ct=C';
At=A';
CA=C*A;
% CtC=Ct*C;
% Ip=speye(p);
% Cv2=zeros(e,1); % initialization needed

%-------------------------------------------------------------------------%
% stuffs for fft-based inversion
%-------------------------------------------------------------------------%
% Circulant matrix to invert via fft
H=(Ct*C)+speye(pp); 

% spectrum of matrix H...ie, the fft of its 1st column
h=fftn(reshape(full(H(:,1)),NSIZE),NSIZE); 
hh = repmat(h,[ones(size(NSIZE)),q]);

%=========================================================================%
% keep track of function value (optional, as it could slow down algorithm)
%=========================================================================%
vec = @(W)W(:);
if MTL
    % last term sums the 2-norm of the rows
    fval = @(W) 1/2 * norm( vec(Y-X*W) )^2 + lam*sum(sqrt(sum(W.^2,2))) + ...
                                           + gam*norm( vec(B.*(CA*W)), 1);
else
    fval = @(W) 1/2 * norm( vec(Y-X*W) )^2 + lam*norm(W(:),1) + ...
                                           + gam*norm( vec(B.*(CA*W)), 1);
end
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
w_old=W;
% disp('go')

% keep track of function value
if funcval, fvalues=zeros(maxiter,1); end;
if exist('wtrue','var'), wdist=zeros(maxiter,1); end;
for k=1:maxiter
%     keyboard
    if funcval,  fvalues(k)=fval(W); end;   
    if exist('wtrue','var'), wdist(k)=norm(W(:)-wtrue(:)); end;
    
    if mod(k,progress)==0 && k~=1
        str='--- %3d out of %d ... Tol=%2.2e (tinner=%4.3fsec, ttotal=%4.3fsec) ---\n';
        fprintf(str,k,maxiter,rel_change,toc(time.inner),toc(time.total))
        time.inner = tic;
    end    
    
    %======================================================================
    % update first variable block: (w,v2)
    %======================================================================
    % update w (if p > n, apply inversion lemma)
%     keyboard
    Q = Xty + rho*(V1-U1) + rho*(At*(V3-U3));
    if p > n
        W=Q/(2*rho) - 1/(2*rho)^2*(K*(X*Q));
    else
        W = (XtX + rho*Ip)\Q;
    end
    
    % update v2
%     keyboard
    V2=C*V3-U2;
    V2(b,:)=soft(V2(b,:),gam/rho);  
%     V2 = bsoft(C*V3 - U2, gam/rho);
%     keyboard

    %======================================================================
    % update second variable block: (v1,v3)
    %======================================================================
    % update v1 
    if MTL
        V1 = vsoft(W+U1,lam/rho);
    else
        V1 = soft(W+U1,lam/rho);
    end   
    
    % update v3 (use fft)
    TMP=Ct*(V2+U2)+A*W+U3;
%     for ii=1:q
% %         tmp=(Ct*(V2(:,ii)+U2(:,ii)))+(A*W(:,ii)+U3(:,ii));
% %         tmp= reshape(tmp, NSIZE);
%         tmp= reshape(TMP(:,ii), NSIZE);
%         tmp=ifftn( fftn(tmp,NSIZE)./h);
%         V3(:,ii)=tmp(:);
%     end
    %---------------------------------------------------------------------%
    % better loopless version using hh = repmat(h,[1,1,20]);
    %---------------------------------------------------------------------%
    TMP = reshape(TMP, [NSIZE,q]);%keyboard
    TMP = ifftn(fftn(TMP)./hh); %keyboard
    V3 = reshape(TMP, [pp,q]);        
    
    %======================================================================
    % dual updates
    %======================================================================
    U1=U1+(W-V1);
    U2=U2+(V2-C*V3);
    U3=U3+(A*W-V3);

    %======================================================================
    % Check termination criteria
    %======================================================================
    %%% relative change in primal variable norm %%%
    rel_change=norm(W(:)-w_old(:))/norm(w_old(:));
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
    w_old=W;
end

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
