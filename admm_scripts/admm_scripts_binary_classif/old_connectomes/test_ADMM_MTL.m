%% test_ADMM_MTL.m
%-------------------------------------------------------------------------%
% (02/14/2014)
% - MTL version
%-------------------------------------------------------------------------%
%%
clear all; 
close all;
drawnow
%% load data
dataOption.datapath=[get_rootdir,'/data_local/ADHD/matfiles/', ...
                                 'designMatrix_adhd_censored.mat'];
% dataOption.datapath=[get_rootdir,'/data_local/autism/matfiles/', ...
%                                  'designMatrix_autism_ABIDE_censored.mat'];
dataOption.dataVars={'X','y','data_info'};
load(dataOption.datapath,dataOption.dataVars{:})
% data_info
[n,p]=size(X);

% optional: save timing result
% saveTime = true;
saveTime = false;
if saveTime
    timingOutputPath = [fileparts(mfilename('fullpath')), '/timing_results/',mfilename];
    timingOutputPath = [timingOutputPath];
end
%% organize site info
% number of sites/task
T = data_info.nSites;
site_mask = cell(T,1);
nlist = zeros(T,1);
Xlist = cell(T,1);
ylist = cell(T,1);

site_info = cell(T+1,4); % <- not really necessary...just for my own convenience
site_info{1,1}='n';
site_info{1,2}='n_DS';
site_info{1,3}='n_HC';
site_info{1,4}='---site---';
for t=1:T
    site_mask{t}=data_info.siteCode == t;
    Xlist{t}=X( site_mask{t}, :);
    ylist{t}=y(site_mask{t});
    nlist(t)=length(ylist{t});
    
    site_info{t+1,1}=length(ylist{t});
    site_info{t+1,2}=sum(ylist{t}==+1);
    site_info{t+1,3}=sum(ylist{t}==-1);
    site_info{t+1,4}=data_info.siteList{t};
end
% site_info
% Xlist,ylist,options,Klist
% return
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
% penalty = 'enet'; % {'enet','gnet','flas','isoTV'}
% penalty = 'gnet';
% penalty = 'flas';
penalty = 'isoTV';
options.lambda=2^-3; % L1 penalty weight
options.gamma =2^-3; % second penalty weight

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
options.termin.progress = 5;   % <- display "progress" (every k iterations...set to inf to disable)
options.termin.silence = true; % <- display termination condition

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
%-------------------------------------------------------------------------%
% optional: save timing result
%-------------------------------------------------------------------------%
if saveTime
    % profiler output path
    timingOutputPath = [timingOutputPath, '_',penalty,'_faster']
    % start profiler
    profile on
end

switch penalty
    case 'enet'
        output=tak_admm_mtl_elasticnet(Xlist,ylist,options);
    case 'gnet'
        output=tak_admm_mtl_graphnet(Xlist,ylist,options);
    case 'flas'
        output=tak_admm_mtl_fusedlasso(Xlist,ylist,options);
    case 'isoTV'
        output=tak_admm_mtl_isotropicTV(Xlist,ylist,options);
end

%-------------------------------------------------------------------------%
% optional: save timing result
%-------------------------------------------------------------------------%
if saveTime
    profile viewer
    profsave(profile('info'),timingOutputPath)
end
%%
% purge
% W=output.W;
W=output.V2;

% number of iteration it took to converge
k=output.k;

% nnz_w=nnz(w)
nnz_W=sum(W~=0)

% training accuracy (not veyr meaningful, but sanity check)
ypred=[];
ytrue=[];
for t=1:T
    ypred=[ypred; SIGN(Xlist{t}*W(:,t))];
    ytrue=[ytrue; ylist{t}];
end
tak_binary_classification_summary(ypred,ytrue)
%% look at the network structure of the support of the parameter
%-------------------------------------------------------------------------%
% get yeo-label info
%-------------------------------------------------------------------------%
load yeo_info347_dilated5mm.mat yeoLabels roiMNI roiLabel
% circularly shift 1 indices (so "unlabeled" is at the final label index)
roiLabel=roiLabel-1;
roiLabel(roiLabel==-1)=12;
[idxsrt,labelCount] = tak_get_yeo_sort(roiLabel);

wlist=cell(T,1);
wmatList=cell(T,1);
wmatSortList=cell(T,1);
for t=1:T
    wlist{t}        = W(:,t);
    wmatList{t}     = tak_dvecinv(W(:,t),0);
    wmatSortList{t} = wmatList{t}(idxsrt,idxsrt);
end

% visualization of the result
purge
cbarOption={'fontsize',22','fontweight','b','ytick',[-.66,0,.66],...
    'YTickLabel',{' <0',' =0',' >0'},'TickLength',[0 0]};
textOption1={'fontweight','b','fontsize',11};
lineOption = {'color','k','linewidth',0.5};
lwidth_deg=2.5;

figure(1),imexp
figure(2),imexp
% figure(3),imexp
for t=1:T
    %---------------------------------------------------------------------%
    % for some damn reason, axis ij needed above...fuck matlab...
    %---------------------------------------------------------------------%
    figure(1),    subplot(2,4,t),hold on
    imcov(wmatSortList{t}),axis off ij,% colorbar(cbarOption{:})
    colorbar off
    tmp=max(abs(caxis));
    caxis([-tmp,tmp]/5)
%     colorbar('location','northoutside')
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
    
    figure(2),    subplot(2,4,t),hold on
    imtriag(wmatSortList{t}),axis off ij,% colorbar(cbarOption{:})
    tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)
end
drawnow