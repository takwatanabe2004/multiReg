function output=tak_admm_fusedlasso_MM(X,y,options,K)
% output=tak_admm_fusedlasso_MM(X,y,options,K)
% (05/13/2014)
%=========================================================================%
% - Combine FC and sMRI features...with spatial penalties
%=========================================================================%
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
    K=tak_admm_inv_lemma(X,1/2);
end

%==========================================================================
% information needed for data augmentation and fft tricks
%==========================================================================
%-------------------------------------------------------------------------%
% the first part: (FC part)
%-------------------------------------------------------------------------%
NSIZE1=options.misc.NSIZE1; % size of the 6d coordinate system
A1=options.misc.A1; % augmentation matrix
At1=A1';
pp1=size(A1,1);

if isfield(options.misc, 'C')
    C1=options.misc.C1; % difference matrix
else
    C1=tak_diffmat(options.misc.NSIZE1,1); % <- C'*C has circulant structure!
end
e1=size(C1,1);
Ct1=C1';

H=(Ct1*C1)+speye(pp1); % Circulant matrix to invert via fft
h1=fftn(reshape(full(H(:,1)),NSIZE1),NSIZE1); % spectrum of matrix H...ie, the fft of its 1st column

b1=logical(options.misc.b1); % masking vector
%-------------------------------------------------------------------------%
% the second part: (sMRI part)
%-------------------------------------------------------------------------%
NSIZE2=options.misc.NSIZE2; % size of the 3d coordinate system
A2=options.misc.A2; % augmentation matrix
At2=A2';
pp2=size(A2,1);

if isfield(options.misc, 'C')
    C2=options.misc.C2; % difference matrix
else
    C2=tak_diffmat(options.misc.NSIZE2,1); % <- C'*C has circulant structure!
end
e2=size(C2,1);
Ct2=C2';

H=(Ct2*C2)+speye(pp2); % Circulant matrix to invert via fft
h2=fftn(reshape(full(H(:,1)),NSIZE2),NSIZE2); % spectrum of matrix H...ie, the fft of its 1st column

b2=logical(options.misc.b2); % masking vector

%-------------------------------------------------------------------------%
% combined part
%-------------------------------------------------------------------------%
A = blkdiag(A1,A2);
%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
[n,p]=size(X);
p1=options.misc.p1; % size of FC
p2=options.misc.p2; % size of sMRI

% primal variable
w =zeros(p,1); 
v1=zeros(n,1);
v2=zeros(p,1);
v3a=zeros(e1,1);
v3b=zeros(e2,1);
v4a=zeros(pp1,1);
v4b=zeros(pp2,1);
v4=zeros(pp1+pp2,1);

% dual variables
u1=zeros(n,1);
u2=zeros(p,1);
u3a=zeros(e1,1);
u3b=zeros(e2,1);
u4a=zeros(pp1,1);
u4b=zeros(pp2,1);
u4=zeros(pp1+pp2,1);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
YX=diag(y)*X;
YXt=YX';
Cv4a=zeros(e1,1); % initialization needed
Cv4b=zeros(e2,1); % initialization needed

% for simplifying A'*(v4-u4) during w update
% mask1a=logical(sum(A1,2));
% mask1b=logical(sum(A2,2));
mask=logical(sum(A,2));

% the diagonal component for updating v3
dv3a=gamma*b1+rho*ones(length(b1),1); 
dv3b=gamma*b2+rho*ones(length(b2),1); 
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
    % update first variable block: (w,v3)
    %======================================================================
    % update w
    q=(YXt*(v1-u1))+(v2-u2)+(v4(mask)-u4(mask));
    w=q/2 - (1/2^2)*(K*(X*q));

    % update v3 = [v3a;v3b]
    v3a=Cv4a-u3a;
    v3b=Cv4b-u3b;
    v3a(b1)=soft(v3a(b1),gamma/rho);  
    v3b(b2)=soft(v3b(b2),gamma/rho);  
%     v3=[v3a;v3b];
    
    %%% precompute terms used more than once %%%
    YXw=YX*w;
%     Aw=A*w;
    Aw1=A1*w(1:p1);
    Aw2=A2*w(1+p1:p1+p2);

    %======================================================================
    % update second variable block: (v1,v2,v4)
    %======================================================================
    % update v1
    v1 = loss.prox(YXw+u1,1/(n*rho));
        
    % update v2
    v2=soft(w+u2,lambda/rho);
    
    %*********************************************************************%
    % update v4=[v4a; v4b]: use fft!
    %*********************************************************************%
    %-----------------------------------%
    % update first block
    %-----------------------------------%
    tmp=(Ct1*(v3a+u3a))+(Aw1+u4a);
    tmp= reshape(tmp, NSIZE1);
    v4a=ifftn( fftn(tmp,NSIZE1)./h1, NSIZE1);
    v4a=v4a(:);
    Cv4a=C1*v4a;
    
    %-----------------------------------%
    % update second block
    %-----------------------------------%
    tmp=(Ct2*(v3b+u3b))+(Aw2+u4b);
    tmp= reshape(tmp, NSIZE2);
    v4b=ifftn( fftn(tmp,NSIZE2)./h2, NSIZE2);
    v4b=v4b(:);
    Cv4b=C2*v4b;
    
    v4=[v4a;v4b];
    
    %======================================================================
    % dual updates
    %======================================================================    
    u1=u1+(YXw-v1);
    u2=u2+(w-v2);
    u3a=u3a+(v3a-Cv4a);
    u3b=u3b+(v3b-Cv4b);
    u4a=u4a+(Aw1-v4a);
    u4b=u4b+(Aw2-v4b);

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
output.v3a=v3a;
output.v3b=v3b;
output.v4a=v4a;
output.v4b=v4b;

% dual variables
output.u1=u1;
output.u2=u2;
output.u3a=u3a;
output.u3b=u3b;
output.u4a=u4a;
output.u4b=u4b;

% time it took for the algorithm to converge
output.time=time.total;

% number of iteration it took to converge
output.k=k;

% final relative change in the primal variable
output.rel_change=rel_change;

% the "track-record" of the relative change in primal variable
output.rel_changevec=rel_changevec(1:k);