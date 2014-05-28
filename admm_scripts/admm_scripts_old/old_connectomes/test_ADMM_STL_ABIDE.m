%% test_ADMM_STL_ABIDE.m
%-------------------------------------------------------------------------%
% (02/14/2014)
% - run single-task admm on the ABIDE dataset
%-------------------------------------------------------------------------%
% (02/13/2014)
% - now use the new circmat convention (wrap around effect at the 1st row of the diffmat)
%-------------------------------------------------------------------------%
% (02/03/2014)
% - empirical risk scaled by 1/n
% - use relative change primal variable as termination condition
%%
clear all; 
close all;
drawnow
%% load data
dataOption.site='ABIDE'; % {'ABIDE','caltech','kki','leuven','maxmun','nyu,...
%                             'ohsu','sdsu','trinity','ucla','um','usm','yale'}
dataOption.datapath=[get_rootdir,'/data_local/autism_info/matfiles/',...
                    'designMatrix_autism_',dataOption.site,'_censored.mat'];
dataOption.dataVars={'X','y'};
load(dataOption.datapath,dataOption.dataVars{:})
num_Nans=sum(isnan(X(:)))
[n,p]=size(X);
nDS=sum(y==+1)
nHC=sum(y==-1)
% return
% reorder data to cluster together indices for the -1 and +1 group
[y,idx]=sort(y);
X=X(idx,:);
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
penalty = 'enet'; % {'enet','gnet','flas','isoTV'}
% penalty = 'gnet';
options.lambda=2^-9; % L1 penalty weight
options.gamma =2^-3; % second penalty weight

% penalty = 'flas';
penalty = 'isoTV';
options.lambda=2^-8.5; % L1 penalty weight
options.gamma =2^-10.5; % second penalty weight

% lambda=2^-9;
% rho = 0.5;
% options.lambda=lambda*rho; % L1 penalty weight
% options.gamma =lambda*(1-rho); % second penalty weight

%==========================================================================
% augmented lagrangian parameters
%==========================================================================
options.rho=1;

%==========================================================================
% termination criterion
%==========================================================================
options.termin.maxiter = 400;   % <- maximum number of iterations
options.termin.tol = 5e-3;      % <- relative change in the primal variable
options.termin.progress = 25;   % <- display "progress" (every k iterations...set to inf to disable)
options.termin.silence = false; % <- display termination condition

%==========================================================================
% information needed for data augmentation and fft tricks
%==========================================================================
if ~strcmpi(penalty,'enet')
    load graph_info347 adjmat coord
    load augmat_mask347newcirc A b

    options.misc.NSIZE=[coord.NSIZE,coord.NSIZE];
    options.misc.A=A; % <- augmentation matrix
    options.misc.b=b; % <- masking vector

    % preassigning this matrix C is optional: 
    % the admm script will compute it internally if the field doesn't exist
    % (helpful when doing CV, since creating this diffmat has minor overhead)
    C=tak_diffmat_newcirc(options.misc.NSIZE,1); % <- C'*C has circulant structure!
    options.misc.C=C;
end
%% run ADMM
switch penalty
    case 'enet'
        output=tak_admm_stl_elasticnet(X,y,options);
    case 'gnet'
        output=tak_admm_stl_graphnet(X,y,options);
    case 'flas'
        output=tak_admm_stl_fusedlasso(X,y,options);
    case 'isoTV'
        output=tak_admm_stl_isotropicTV(X,y,options);
end
%%
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

% figure,imexp
% % primal variables and residuals
% subplot(261),tplot(log10(output.rel_changevec))
% subplot(268),tplot(y)
% subplot(269),tplot(ypred)
% subplot(264),tplot(w)
% subplot(265),tplot(v1),% ylim([-.6 .6])
% subplot(266),tplot(v2),% ylim([-.6 .6])
% subplot(2,6,11),tplot(YX*w-v1),title('resid1 (YXw-v1)')
% subplot(2,6,12),tplot(w-v2),title('resid2 (w-v2)')
% % subplot(2,6,11)
% % subplot(2,6,12)
% %%%
% % (y==(ypred))'
% % tak_binary_classification_summary(ypred,y)
% % %%
% % (yts==-sign(Xts*w))'
%% look at the network structure of the support of the parameter
load yeo_info347_dilated5mm.mat yeoLabels roiMNI roiLabel
[idx_sort,labelCount] = tak_get_yeo_sort(roiLabel);

v2mat=tak_dvecinv(v2)-eye(347);
v2matsort=v2mat(idx_sort,idx_sort);

% figure,imexp
% subplot(121),imcov(v2mat)
% subplot(122),imcov(v2matsort)
% tak_box_covGroups(gcf, labelCount, yeoLabels)
% figure,imexp
% subplot(121),imtriag(v2mat)
% subplot(122),imtriag(v2matsort)
% tak_box_covGroups(gcf, labelCount, yeoLabels)
%% (12/24/2013)show the figure i created for the paper
% circularly shift 1 indices
roiLabel=roiLabel-1;
roiLabel(roiLabel==-1)=12;
% return

[idxsrt,labelCount] = tak_get_yeo_sort(roiLabel);
v2matsort2=v2mat(idxsrt,idxsrt);
%%
cbarOption={'fontsize',22','fontweight','b','ytick',[-.66,0,.66],...
    'YTickLabel',{' <0',' =0',' >0'},'TickLength',[0 0]};
lwidth_deg=2.5;
figure,imexpl,imtriag(v2mat),axis off,colorbar(cbarOption{:})
figure,imexpl
imtriag(v2matsort2),axis off
% tak_box_covGroups(gcf, labelCount, yeoLabels)
deg=sum(v2matsort2~=0);
hold on
plot(347*ones(1,347)-deg,'k','linewidth',lwidth_deg)
textOption1={'fontweight','b','fontsize',12};
lineOption = {'color','k','linewidth',0.5};
colorbar(cbarOption{:})
tak_local_linegroups3(gcf,labelCount,textOption1,yeoLabels,lineOption),hold on

figure,imexpl
imcov(v2matsort2),,axis off
tmp=max(abs(caxis));
caxis([-tmp,tmp]/5)
% tak_box_covGroups(gcf, labelCount, yeoLabels)
deg=sum(v2matsort2~=0);
hold on
plot(347*ones(1,347)-deg,'k','linewidth',lwidth_deg)
tak_local_linegroups3(gcf,labelCount,textOption1,yeoLabels,lineOption),hold on
% save enet_weight_lam_minus6_gam_minus1 v2matsort2 labelCount textOption1 yeoLabels lineOption options
%%
% purge
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