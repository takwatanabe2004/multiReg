function output=tak_admm_stl_isotropicTV(Xlist,ylist,options,Klist)
% output=tak_admm_mtl_isotropicTV(Xlist,ylist,options,Klist)
% (02/22/2014)
%-------------------------------------------------------------------------%
% Xlist: (Tx1) cell array containing T design matrices (one for each site)
%  yist: (Tx1) cell array containing T label vector    (one for each site)
% Klist: (Tx1) cell array containing T (nt x p) matrices for the inversion
%              lemma (one for each site)
%-------------------------------------------------------------------------%
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
% information needed for data augmentation and fft tricks
%==========================================================================
NSIZE=options.misc.NSIZE; % size of the 6d coordinate system
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

%=========================================================================%
% multitask info
%=========================================================================%
T = length(ylist); % number of sites/tasks
YXlist =cell(T,1);  
YXtlist=cell(T,1);  
nlist = zeros(T,1); % number of subjects in each site
for t=1:T
    YXlist{t}=bsxfun(@times,Xlist{t},ylist{t}); % bsxfun faster than diag(y)*X;
    YXtlist{t}=YXlist{t}';
    nlist(t)=length(ylist{t});
end
p=size(Xlist{1},2);

%-------------------------------------------------------------------------%
% Matrix K for inversion lemma
%-------------------------------------------------------------------------%
if nargin < 4 % compute matrix "K" needed for inversion lemma
    Klist = cell(T,1);
    for t=1:T
        Klist{t}=tak_admm_inv_lemma(Xlist{t},1/2);
    end
end
%% initialize variables, function handles, and terms used through admm steps
%==========================================================================
% initialize variables
%==========================================================================
% primal variable
W =zeros(p,T); 
V1=cell(T,1); % <- use cell since #subjects differ by each site
V2=zeros(p,T);
V3=zeros(e,T);
V4=zeros(pp,T);

% dual variables
U1=cell(T,1); % <- use cell since #subjects differ by each site
U2=zeros(p,T);
U3=zeros(e,T);
U4=zeros(pp,T);

for t=1:T
    V1{t}=zeros(nlist(t),1);
    U1{t}=zeros(nlist(t),1);
end

%==========================================================================
% function handles
%==========================================================================
soft=@(t,tau) sign(t).*max(0,abs(t)-tau); % soft-thresholder

%=========================================================================%
% initialize terms used throughout the admm iterations
%=========================================================================%
YXw_list = cell(T,1);
CV4=zeros(e,T); % initialization needed

% for simplifying A'*(v4-u4) during w update
mask1=logical(sum(A,2));

%-------------------------------------------------------------------------%
% create index cell array of the coordinate system
% - needed for the prox operator unique to the isotropic tv norm (v3 update)
%-------------------------------------------------------------------------%
DIM = length(NSIZE);
idxCell=cell(DIM,1);
for i=1:DIM
    idxCell{i}=1+(i-1)*pp:i*pp;
end
%% begin admm iteration
time.total=tic;
time.inner=tic;

rel_changevec=zeros(maxiter,1);
W_old=W;
% disp('go')
for k=1:maxiter
    if mod(k,progress)==0 && k~=1        
        str='--- %3d out of %d ... Tol=%2.2e (tinner=%4.3fsec, ttotal=%4.3fsec) ---\n';
        fprintf(str,k,maxiter,rel_change,toc(time.inner),toc(time.total))
        time.inner = tic;
    end
    
    %======================================================================
    % update first variable block: (W,V3)
    %======================================================================
    for t=1:T
        %-----------------------------------------------------------------%
        % update W
        %-----------------------------------------------------------------%
        q =   YXtlist{t} * (V1{t}-U1{t}) ...
            + V2(:,t)-U2(:,t)  +  V4(mask1,t)-U4(mask1,t);
        W(:,t) = q/2 - (1/4) * (Klist{t}*(Xlist{t}*q));

        %-----------------------------------------------------------------%
        % precompute terms used more than once
        %-----------------------------------------------------------------%
        YXw_list{t}=YXlist{t}*W(:,t);
        
        %-----------------------------------------------------------------%
        % update V3
        %-----------------------------------------------------------------%
%         V3(:,t) = tak_isoTV_prox(CV4(:,t)-U3(:,t),b,gamma/rho,idxCell);
    end
    %---------------------------------------------------------------------%
    % update V3
    % - this approach is faster than looping through the function 
    %   tak_isoTV_prox.m for each task 
    %   (obviates the need to loop through the tasks)
    % - note: script only considers 6-D functional connectome (no 4-D)
    %---------------------------------------------------------------------%
    %*****************************************************************%
    % Signal which the prox-operator of the istropic TV will be applied
    %*****************************************************************%
    signal=CV4-U3;

    % apply masking to avoid artifacts from data-augmentation & circulancy
    sigMasked=bsxfun(@times,b,signal);

    % the gradient of the (masked) signal in the 6 coordinate directions
    DX1=sigMasked(idxCell{1},:);
    DY1=sigMasked(idxCell{2},:);
    DZ1=sigMasked(idxCell{3},:);
    DX2=sigMasked(idxCell{4},:);
    DY2=sigMasked(idxCell{5},:);
    DZ2=sigMasked(idxCell{6},:);

    %*********************************************************************%
    % euclidean norm of the gradients 
    % (the sum of this is the isotropic TV of the signal)
    %*********************************************************************%
    sigNormMasked = sqrt(DX1.^2 + DY1.^2 + DZ1.^2 + DX2.^2 + DY2.^2 + DZ2.^2);

    %*********************************************************************%
    % Compute shrinkage factor with the convention (0/0)*0 = 0
    %*********************************************************************%
    tmp=(gamma/rho)./sigNormMasked;
    tmp(isnan(tmp))=0;
    shrink_f=max(1-tmp,0); % shrinkage factor

    %*********************************************************************%
    % replicate shrinkage factor for the 6 coordinate directions
    % (faster than "repmat" or looping, plus "bsxfun" does not apply here)
    %*********************************************************************%
    shrink=[shrink_f;shrink_f;shrink_f;shrink_f;shrink_f;shrink_f];

    %*********************************************************************%
    % apply shrinkage
    %*********************************************************************%
    shrinked_signal=shrink.*signal;

    %*********************************************************************%
    % The result of the prox operator has two components:
    % (1) part where b==0 is left unchanged (returns the trivial solution)
    % (2) part where b==1 receives the shrinked signal component
    %*********************************************************************%
    V3=signal;
    V3(b,:)=shrinked_signal(b,:);

    %---------------------------------------------------------------------%
    % precompute terms used more than once 
    %---------------------------------------------------------------------%
    AW=A*W;

    %======================================================================
    % update second variable block: (V1,V2,V4)
    %======================================================================
    % temp array for fft inside the loop
    TMP = (Ct*(V3+U3)) + (AW+U4);
    for t=1:T
        %-----------------------------------------------------------------%
        % update V1
        %-----------------------------------------------------------------%
        V1{t} = loss.prox(YXw_list{t} + U1{t}, 1/(nlist(t)*rho));

        %-----------------------------------------------------------------%
        % update V4
        %-----------------------------------------------------------------%
        tmp = reshape(TMP(:,t), NSIZE);
        tmp = ifftn( fftn(tmp,NSIZE)./h);
        V4(:,t) = real(tmp(:));        
    end
    CV4=C*V4;
    
    %---------------------------------------------------------------------%
    % update V2
    %---------------------------------------------------------------------%
%     V2=tsoftvec(W+U2,lambda/rho);
    V2=tak_soft(W+U2,lambda/rho); % <- single task case
        
    %======================================================================
    % dual updates
    %======================================================================  
    for t=1:T
        U1{t}   = U1{t}   + (YXw_list{t}-V1{t});
    end
    U2 = U2 + (W-V2);
    U3 = U3 + (V3-CV4);
    U4 = U4 + (AW-V4);
    
    %======================================================================
    % Check termination criteria
    %======================================================================
    %%% relative change in primal variable norm %%%
    rel_change=norm(W-W_old,'fro')/norm(W_old,'fro');
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
    W_old=W;
end
time.total=toc(time.total);
%% organize output
% primal variables
output.W=W;
output.V1=V1;
output.V2=V2;
output.V3=V3;
output.V4=V4;

% dual variables
output.U1=U1;
output.U2=U2;
output.U3=U3;
output.U4=U4;

% time it took for the algorithm to converge
output.time=time.total;

% number of iteration it took to converge
output.k=k;

% final relative change in the primal variable
output.rel_change=rel_change;

% the "track-record" of the relative change in primal variable
output.rel_changevec=rel_changevec(1:k);