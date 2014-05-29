% sim_run_admm.m (01/07/2014)
%==========================================================================
%--------------------------------------------------------------------------
%%
clear
% purge

snr=0.4;
nn_struct = 4;

if nn_struct == 4
    dataPath=[fileparts(mfilename('fullpath')),'/sim_dataset_snr',num2str(snr),'.mat'];
elseif nn_struct==8
    dataPath=[fileparts(mfilename('fullpath')),'/sim_8nn_dataset_snr',num2str(snr),'.mat'];
else
    error('4nn or 8nn are the only options')
end
dataVars={'sim_info','Xtr','ytr','Xts','yts'};
load(dataPath,dataVars{:})
%% select number of training samples (ceiling: sim_info.ntr)
ntr = 200;
% ntr = sim_opt.ntr;
idx_plus=1:ntr/2;
idx_minu=sim_info.ntr/2+1:sim_info.ntr/2+ntr/2;
idx=[idx_plus,idx_minu];
ytr=ytr(idx);
Xtr=Xtr(idx,:);
% [sum(ytr==+1),sum(ytr==-1)]
% return
%% set penalty parameters
% penalty = 'enet'; % {'enet','gnet','flas'}
% penalty = 'gnet';
% penalty = 'flas';
penalty = 'isoTV';

% options.lambda=2^-6;  % L1 penalty weight
% options.gamma =2^-1; % fused lasso penalty weight

options.lambda=2^-8;  % L1 penalty weight
options.gamma =2^-7.25; % fused lasso penalty weight

%% set algorithm options (this block doesn't need to be touched)
% loss function
options.loss='hinge1';

% augmented lagrangian parameters
options.rho=1;

% termination criterion
options.termin.maxiter = 400;   % <- maximum number of iterations
options.termin.tol = 4e-3;      % <- relative change in the primal variable
options.termin.progress = 500;   % <- display "progress" (every k iterations...set to inf to disable)
options.termin.silence = false; % <- display termination condition

% information needed for data augmentation and fft tricks
load graph_info347_2d adjmat coord
load Amat_Bmat347_2d.mat A b

options.misc.NSIZE=[coord.NSIZE,coord.NSIZE];
options.misc.A=A; % <- augmentation matrix
options.misc.b=b; % <- masking vector
%% run ADMM
switch penalty
    case 'enet'
        output=tak_admm_elasticnet(Xtr,ytr,options);
    case 'gnet'
        output=tak_admm_graphnet(Xtr,ytr,options);
    case 'flas'
        output=tak_admm_fusedlasso(Xtr,ytr,options);
    case 'isoTV'
        output=tak_admm_isotropicTV(Xtr,ytr,options);
end
w=output.v2;

training_result=tak_binary_classification_summary(SIGN(Xtr*w),ytr)
ypr=SIGN(Xts*w);
tak_binary_classification_summary(ypr,yts)
% return
anom_nodes_mask = sim_info.anom_nodes.mask;
% imcovvl(tak_dvecinv(w,0))

figure,imexpb
subplot(131),imcov(tak_dvecinv(w,0)), colorbar off
subplot(132),imcov(abs(tak_dvecinv(w,0))),colorbar off
subplot(133),imcov(anom_nodes_mask),colorbar off

% figure,imexpb
% subplot(132),imedge(tak_dvecinv(w,0))
% subplot(133),imedge(anom_nodes_mask)
% subplot(131),imcov(tak_dvecinv(w,0)),colormap('jet'),
% % colorbar location NorthOutside
% % colorbar off