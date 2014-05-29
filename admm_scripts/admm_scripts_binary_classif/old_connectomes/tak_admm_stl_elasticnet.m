function output=tak_admm_stl_elasticnet(Xlist,ylist,options,Klist)
% output=tak_admm_stl_elasticnet(Xlist,ylist,options,Klist)
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
        Klist{t}=tak_admm_inv_lemma(Xlist{t},rho/(rho+gamma));
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

% dual variables
U1=cell(T,1); % <- use cell since #subjects differ by each site
U2=zeros(p,T);

for t=1:T
    V1{t}=zeros(nlist(t),1);
    U1{t}=zeros(nlist(t),1);
end

%=========================================================================%
% initialize terms used throughout the admm iterations
%=========================================================================%
YXw_list = cell(T,1);
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
    % update first variable block: (W)
    %======================================================================
    % update W
    for t=1:T
        q =   rho*(  YXtlist{t} * (V1{t}-U1{t})  ) ...
            +                 rho*(   V2(:,t)-U2(:,t)  );
        W(:,t) = q/(rho+gamma) ...
                 - (rho/(rho+gamma)^2) * (Klist{t}*(Xlist{t}*q));
    end
      
    %%% precompute terms used more than once %%%
    for t=1:T
        YXw_list{t}=YXlist{t}*W(:,t);
    end

    %======================================================================
    % update second variable block: (V1,V2)
    %======================================================================
    % update V1
    for t=1:T
        V1{t} = loss.prox(YXw_list{t} + U1{t}, 1/(nlist(t)*rho));
    end
        
    % update v2
%     V2=tsoftvec(W+U2,lambda/rho);
    V2=tak_soft(W+U2,lambda/rho); % <- single task case
    
    %======================================================================
    % dual updates
    %======================================================================
    for t=1:T
        U1{t}=U1{t} + (YXw_list{t}-V1{t});
        U2(:,t)=U2(:,t) + (W(:,t)-V2(:,t));
    end

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

% dual variables
output.U1=U1;
output.U2=U2;

% time it took for the algorithm to converge
output.time=time.total;

% number of iteration it took to converge
output.k=k;

% final relative change in the primal variable
output.rel_change=rel_change;

% the "track-record" of the relative change in primal variable
output.rel_changevec=rel_changevec(1:k);