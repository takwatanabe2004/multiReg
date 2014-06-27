% sim_run_admm_gnet.m (01/07/2014)
%==========================================================================
%--------------------------------------------------------------------------
%%
clear
purge

snr=0.5;
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
%% set penalty parameters
% options.lambda=0;  % L1 penalty weight
% options.gamma =0; % graphnet penalty weight

options.lambda=2^-9;  % L1 penalty weight
options.gamma =2^-10; % graphnet penalty weight
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
load augmat_mask347ver2_2d A b

options.misc.NSIZE=[coord.NSIZE,coord.NSIZE];
options.misc.A=A; % <- augmentation matrix
options.misc.b=b; % <- masking vector
%% run ADMM
output=tak_admm_graphnet(Xtr,ytr,options);

w=output.v2;

ypr=SIGN(Xts*w);
tak_binary_classification_summary(ypr,yts)

anom_nodes_mask = sim_info.anom_nodes.mask;
imcovvl(tak_dvecinv(w,0))

figure,imexpb
subplot(121),imedge(tak_dvecinv(w,0))
subplot(122),imedge(anom_nodes_mask)