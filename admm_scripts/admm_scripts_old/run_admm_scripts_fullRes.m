%% run_admm_scripts_fullRes.m
% (05/10/2014)
%=========================================================================%
% - run admm scripts on full data (no CV procedure)
%-------------------------------------------------------------------------%
% - script from may10_run_admm_flas.m
%=========================================================================%
%%
clear
% purge
%% load data
load sMRI_design_censor_fullRes.mat X y

[n,p]=size(X);

% reorder data to cluster together indices for the -1 and +1 group
% [y,idx]=sort(y);
% X=X(idx,:);

% scale to 0 mean, unit variance?
% flag_scale = true;
flag_scale = false;
if flag_scale
    X=zscore(X);
end

%==========================================================================
% load info for the indexing for the 10-fold-CV partitioning
%==========================================================================
% load Overall Overall
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
% % penalty = 'gnet';
% options.lambda=2^-6; % L1 penalty weight
% options.gamma =2^2; % second penalty weight

% penalty = 'flas';
penalty = 'isoTV';
options.lambda=2^-9.0; % L1 penalty weight
options.gamma =2^-8.5; % second penalty weight

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
    
    % augmentation matrix
    load([get_rootdir,'/data_local/A_matrix_cropped_fullRes.mat'],'A','nx','ny','nz')
    
    % mask vector
    load([get_rootdir,'/data_local/B_matrix_cropped_fullRes.mat'],'b')
    %=====================================================================%
%     load graph_info347 adjmat coord
%     load augmat_mask347 A b

    options.misc.NSIZE=[nx,ny,nz];
    options.misc.A=A; % <- augmentation matrix
    options.misc.b=b; % <- masking vector

    % preassigning this matrix C is optional: 
    % the admm script will compute it internally if the field doesn't exist
    % (helpful when doing CV, since creating this diffmat has minor overhead)
    C=tak_diffmat(options.misc.NSIZE,1); % <- C'*C has circulant structure!
    options.misc.C=C;
end
%% run ADMM
switch penalty
    case 'enet'
        output=tak_admm_elasticnet(X,y,options);
    case 'gnet'
        output=tak_admm_graphnet(X,y,options);
    case 'flas'
        output=tak_admm_fusedlasso(X,y,options);
    case 'isoTV'
        output=tak_admm_isotropicTV(X,y,options);
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

%-------------------------------------------------------------------------%
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
tak_binary_classification_summary(ypred,y)
%% visualization
% purge
load([get_rootdir, '/data_local/A_matrix_fullFOV_fullRes.mat'], 'A_full')
nii =load_nii([get_rootdir,'/data_local/m0wrp1mprage.nii']);
[nx,ny,nz]=size(nii.img);

%=========================================================================%
% create nii of the weight vector in the "volume space"
%=========================================================================%
wvol = reshape(A_full*output.v2, [nx,ny,nz]);

%-------------------------------------------------------------------------%
% set the regions outside the brain support to have the minimum 
% intensity value (improves visualization)
%-------------------------------------------------------------------------%
bmask=load_nii([get_rootdir,'/data_local/brain_mask.nii']);
bmask.img=logical(bmask.img);
% view_nii(bmask);
wvol(~bmask.img)=min(wvol(:));


%-------------------------------------------------------------------------%
% create nii version
%-------------------------------------------------------------------------%
nii_wvol = nii;
nii_wvol.img=wvol;

%=========================================================================%
% create thresholded mask volume of the weight vector
%=========================================================================%
%-------------------------------------------------------------------------%
% simple threshold by "epsilon"
%-------------------------------------------------------------------------%
% wvol_supp = zeros(size(wvol));
% wvol_supp(abs(wvol)>1e-3)=+1;

%-------------------------------------------------------------------------%
% threshold top "k" weights with largest magnitude
%-------------------------------------------------------------------------%
% purge
k=20000;
k=min(k,nnz(v2)); % <- incase k exceeds nnz of the weight vector
[w_rank,idx_rank]=sort(abs(v2),'descend');
tmp=zeros(size(v2));
tmp(idx_rank(1:k))=+1;
wvol_supp=reshape(A_full*tmp,[nx,ny,nz]);

%-------------------------------------------------------------------------%
% set background to -1
%-------------------------------------------------------------------------%
wvol_supp(~bmask.img)=-1;
% sum(wvol_supp(:)==1) % <- sanity check

%-------------------------------------------------------------------------%
% create nii version
%-------------------------------------------------------------------------%
nii_wvol_supp = nii;
nii_wvol_supp.img=wvol_supp;

%=========================================================================%
% visualize
%=========================================================================%
% view_nii(nii_wvol);
% view_nii(nii_wvol_supp);
vis3d(wvol, 'jet')
vis3d(wvol_supp, 'multi')
%%
% purge
% 
% % sample brain volume
% brain_vol=reshape(A_full*X(1,:)',[nx,ny,nz]);
% nii_brain_vol=nii;
% nii_brain_vol.img=brain_vol;
% view_nii(nii_brain_vol);
% view_nii(nii_wvol);
% 
% % view_nii(load_nii([get_rootdir,'/data_local/m0wrp1mprage.nii']));
% % view_nii(tak_downsample_nii(load_nii([get_rootdir,'/data_local/m0wrp1mprage.nii'])));
% 
% figure
% vol3d('CData',brain_vol)


%%
% % purge
% % figure,vol3d('CData',wvol,'texture','3D')
% % 
% bmask=tak_downsample_nii(load_nii([get_rootdir,'/data_local/brain_mask.nii']));
% bmask.img=logical(bmask.img);
% % % view_nii(bmask);
% % 
% wvol2=wvol;
% wvol2(~bmask.img)=min(wvol(:));
% wvol_tri2=wvol_supp;
% wvol_tri2(~bmask.img)=-2;%min(wvol_tri2(:));
% % nii_wvol2=nii_wvol;
% % nii_wvol2.img=wvol2;
% % view_nii(nii_wvol2);
% % tak_gui_show_slices(wvol2)
% % figure,vol3d('CData',wvol2,'texture','3D')
% % figure,vol3d('CData',wvol2,'texture','3D')
% vis3d(wvol2)
% vis3d(wvol_tri2)