function output=tak_admm_isotropicTV(X,y,options,K)
% output=tak_admm_isotropicTV(X,y,options,K)
% (02/22/2014) - changed name from tak_admm_stl_isotropicTV.m
% (02/13/2014) - use the new circulant matrix form
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
NSIZE=options.misc.NSIZE; % size of the coordinate system
A=options.misc.A; % augmentation matrix
At=A';
pp=size(A,1);

if isfield(options.misc, 'C')
    C=options.misc.C; % difference matrix
else
    C=tak_diffmat_newcirc(options.misc.NSIZE,1); % <- C'*C has circulant structure!
end
e=size(C,1);
Ct=C';

H=(Ct*C)+speye(pp); % Circulant matrix to invert via fft
h=fftn(reshape(full(H(:,1)),NSIZE),NSIZE); % spectrum of matrix H...ie, the fft of its 1st column

b=logical(options.misc.b); % masking vector

%-------------------------------------------------------------------------%
% create index cell array of the coordinate system
% - needed for the prox operator unique to the isotropic tv norm (v3 update)
%-------------------------------------------------------------------------%
DIM = length(NSIZE);
idxCell=cell(DIM,1);
for i=1:DIM
    idxCell{i}=1+(i-1)*pp:i*pp;
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
v3=zeros(e,1);
v4=zeros(pp,1);

% dual variables
u1=zeros(n,1);
u2=zeros(p,1);
u3=zeros(e,1);
u4=zeros(pp,1);

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder

%==========================================================================
% precompute terms used throughout admm
%==========================================================================
YX=diag(y)*X;
YXt=YX';
Cv4=zeros(e,1); % initialization needed

% for simplifying A'*(v4-u4) during w update
mask1=logical(sum(A,2));
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
    q=(YXt*(v1-u1))+(v2-u2)+(v4(mask1)-u4(mask1));
    w=q/2 - (1/2^2)*(K*(X*q));

    %---------------------------------------------------------------------%
    % update v3
    %---------------------------------------------------------------------%
    v3 = tak_isoTV_prox(Cv4-u3,b,gamma/rho,idxCell);
    
    %%% precompute terms used more than once %%%
    YXw=YX*w;
    Aw=A*w;

    %======================================================================
    % update second variable block: (v1,v2,v4)
    %======================================================================
    % update v1
    v1 = loss.prox(YXw+u1,1/(n*rho));
        
    % update v2
    v2=soft(w+u2,lambda/rho);
    
    % update v4: use fft!
    tmp=(Ct*(v3+u3))+(Aw+u4);
    tmp= reshape(tmp, NSIZE);
    v4=ifftn( fftn(tmp,NSIZE)./h, NSIZE);
    v4=v4(:);
    Cv4=C*v4;
    
    %======================================================================
    % dual updates
    %======================================================================    
    u1=u1+(YXw-v1);
    u2=u2+(w-v2);
    u3=u3+(v3-Cv4);
    u4=u4+(Aw-v4);

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
output.v3=v3;
output.v4=v4;

% dual variables
output.u1=u1;
output.u2=u2;
output.u3=u3;
output.u4=u4;

% time it took for the algorithm to converge
output.time=time.total;

% number of iteration it took to converge
output.k=k;

% final relative change in the primal variable
output.rel_change=rel_change;

% the "track-record" of the relative change in primal variable
output.rel_changevec=rel_changevec(1:k);