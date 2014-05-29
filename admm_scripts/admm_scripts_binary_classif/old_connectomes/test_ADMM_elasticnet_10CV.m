%% test_ADMM_elasticnet_10CV.m
% (12/16/2013)
% - run cross validation
% - note that unlike Graphnet/fusedLasso, the matrix K for the inversion
%   lemma is a function of the regularization paramter \gamma, thus
%   this has to be recomputed everytime \gamma changes.
%%
clear all; 
close all;
drawnow
%% load data
load rcorr_design_censor.mat X y
[n,p]=size(X);

%==========================================================================
% load info for the indexing for the 10-fold-CV partitioning
%==========================================================================
load Overall Overall

%==========================================================================
% yeo network coordinate info (for visualization)
%==========================================================================
load yeo_info347_dilated5mm.mat yeoLabels roiMNI roiLabel
[idx_sort,labelCount] = tak_get_yeo_sort(roiLabel);
%% set algorithm options
%==========================================================================
% loss function
%==========================================================================
options.loss='hinge1';
% options.loss='hinge2';
% options.loss='hubhinge';
% options.loss_huber_param=0.2; % <- only needed when using huberized-hinge

%==========================================================================
% set penalty parameters
%==========================================================================
options.lambda=2^-11; % L1 penalty weight
options.gamma =2^-5; % L2 penalty weight

%==========================================================================
% augmented lagrangian parameters
%==========================================================================
options.rho=1;

%==========================================================================
% termination criterion
%==========================================================================
options.termin.maxiter = 500;   % <- maximum number of iterations
options.termin.tol = 4e-3;      % <- relative change in the primal variable
options.termin.progress = inf;   % <- display "progress" (every k iterations...set to inf to disable)
options.termin.silence = false; % <- display termination condition

%==========================================================================
% precompute the matrix K involved in the inversion lemma for updating w
%  - the form of the matrix to invert is: (ilemma1* X'*X + ilemma2*I)^-1
%==========================================================================
% precompute matrix K (optional): 
K=tak_admm_inv_lemma(X,options.rho/(options.rho+options.gamma));
%% begin 10-fold CV
ypredicted = [];
ytrue      = [];
tic_idxCV=tic;
for idxCV = 1:10
    fprintf('***** idxCV = %2d ...%6.3f sec *****\n',idxCV,toc(tic_idxCV))
    %======================================================================
    % 10-fold-CV data partition
    %======================================================================
    mask_ts = Overall.CrossValidFold(:,idxCV);
    mask_tr = ~mask_ts;
    
    Xts = X(mask_ts,:);
    Xtr = X(mask_tr,:);

    yts = y(mask_ts);
    ytr = y(mask_tr);
    
    %======================================================================
    % run ADMM
    %======================================================================
    output=tak_admm_elasticnet(Xtr,ytr,options,K);
    
    %======================================================================
    % analyze admm result
    %======================================================================
%     w=output.w;
    w=output.v2;
    nnz(w)
    % prediction on the testing data
    ypr=SIGN(Xts*w);
    
    tak_binary_classification_summary(ypr,yts)
    ypredicted = [ypredicted;ypr];
    ytrue      = [ytrue; yts];
    
    if idxCV<=3
        %%
        figure,imexp
        subplot(4,2,7),tplot(w)    
        if idxCV==1,AXIS=axis; else axis(AXIS);end;
        subplot(4,2,8),tplot(sort(w))    
        if idxCV==1,AXIS=axis; else axis(AXIS);end;
        w_mat=tak_dvecinv(w,0);
        subplot(4,2,[1,3,5]),imcov(w_mat),
        if idxCV==1,CAXIS1=caxis; else caxis(CAXIS1);end;
        caxis(caxis./[3 3])
        wmat_sort=w_mat(idx_sort,idx_sort);
        subplot(4,2,[2,4,6]),imcov(wmat_sort)
        tak_box_covGroups(gcf, labelCount, yeoLabels)
        if idxCV==1,CAXIS2=caxis; else caxis(CAXIS2);end;
        caxis(caxis./[3 3])
        %%
    end
    drawnow
end % <- idxCV loop

classification_summary=tak_binary_classification_summary(ypredicted,ytrue)
%% run ADMM on full data
output=tak_admm_elasticnet(X,y,options,K);
% purge
w=output.w;
v1=output.v1;
v2=output.v2;

% number of iteration it took to converge
k=output.k;

% nnz_w=nnz(w)
nnz_v2=nnz(v2)

% ypred=sign(X*v2);
ypred=SIGN(X*v2);
YX=diag(y)*X;

figure,imexp
% primal variables and residuals
subplot(261),tplot(log10(output.rel_changevec))
subplot(268),tplot(y)
subplot(269),tplot(ypred)
subplot(264),tplot(w)
subplot(265),tplot(v1),% ylim([-.6 .6])
subplot(266),tplot(v2),% ylim([-.6 .6])
subplot(2,6,11),tplot(YX*w-v1),title('resid1 (YXw-v1)')
subplot(2,6,12),tplot(w-v2),title('resid2 (w-v2)')
% subplot(2,6,11)
% subplot(2,6,12)
%%%
% (y==(ypred))'
% tak_binary_classification_summary(ypred,y)
% %%
% (yts==-sign(Xts*w))'
%% look at the network structure of the support of the parameter
v2mat=tak_dvecinv(v2)-eye(347);
v2matsort=v2mat(idx_sort,idx_sort);

figure,imexp
subplot(121),imcov(v2mat)
subplot(122),imcov(v2matsort)
tak_box_covGroups(gcf, labelCount, yeoLabels)
figure,imexp
subplot(121),imedge(v2mat)
subplot(122),imedge(v2matsort)
tak_box_covGroups(gcf, labelCount, yeoLabels)
%%
% vmax=max(v)
% wmax=max(w)
% tmax=max(t)
% u1max=max(u1)
% u2max=max(u2)
% u3max=max(u3)
%%
% timeStamp=tak_timestamp;
% mFileName=mfilename;
% save elasticnet_admm v2mat options data_option timeStamp mFileName