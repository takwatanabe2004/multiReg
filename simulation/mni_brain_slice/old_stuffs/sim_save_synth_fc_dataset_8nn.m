% sim_save_synth_fc_dataset_8nn (01/07/2014)
% - same as sim_save_synth_fc_dataset.m, but using 8-nearest neighbor for
%   the anomalous nodes
%==========================================================================
%--------------------------------------------------------------------------
%%
clear all
purge
load('Results_tak.mat', 'ConnMean', 'ConnVar')
load('graph_info347_2d.mat', 'coord')
load('sim_anom_node_info_8nn.mat', 'anom_nodes')
%% set options
% set SNR devel (same as cohen's d) http://en.wikipedia.org/wiki/Cohen%27s_d#Cohen.27s_d 
% http://www.3rs-reduction.co.uk/html/6__power_and_sample_size.html
snr=0.4;

% nHC_tr = number of training control samples 
% nDS_tr = number of training patient samples
%   n_tr = number of training samples
nHC_tr = 100;
nDS_tr=nHC_tr;
ntr=nHC_tr+nDS_tr;

% nHC_ts = number of training control samples 
% nDS_ts = number of training patient samples
%   n_ts = number of training samples
nHC_ts = 250;
nDS_ts=nHC_ts;
nts=nHC_ts+nDS_ts;

seedPoint=0;
randn('seed',seedPoint)
%% set save options
fsave=true;
% fsave=false;
% outPath=[fileparts(mfilename('fullpath')),'/sim_8nn_training_ntr',num2str(ntr),...
%          '_nts',num2str(nts),'_snr',num2str(snr),'.mat']
outPath=[fileparts(mfilename('fullpath')),'/sim_8nn_dataset_snr',num2str(snr),'.mat']
     
outVars={'sim_info','Xtr','Xts','ytr','yts','timeStamp','mFileName'};
% return
%% relevant info
ConnMeanVec = tak_dvec(ConnMean);
ConnVarVec  = tak_dvec(ConnVar);
p=length(ConnMeanVec);
%% create dataset set
%==========================================================================
% create -1 (HC) control subjects
%==========================================================================
X_HC_tr = repmat(ConnMeanVec(:)',[nHC_tr,1])+ ...
          repmat(sqrt(ConnVarVec(:))',[nHC_tr,1]).*randn(nHC_tr,p);
X_HC_ts = repmat(ConnMeanVec(:)',[nHC_ts,1])+ ...
          repmat(sqrt(ConnVarVec(:))',[nHC_ts,1]).*randn(nHC_ts,p);

% inverse fisher tx to map to (-1,+1) correlation space
X_HC_tr = tak_fisher_itransformation(X_HC_tr);
X_HC_ts = tak_fisher_itransformation(X_HC_ts);

% Whos,return
%==========================================================================
% create +1 (DS) disease subjects
%==========================================================================
X_DS_tr = repmat(ConnMeanVec(:)',[nDS_tr,1])+ ...
          repmat(sqrt(ConnVarVec(:))',[nDS_tr,1]).*randn(nDS_tr,p);
X_DS_ts = repmat(ConnMeanVec(:)',[nDS_ts,1])+ ...
          repmat(sqrt(ConnVarVec(:))',[nDS_ts,1]).*randn(nDS_ts,p);
      
% add SIGNAL on connections between the two anomalous node clusters
for i=1:length(anom_nodes.idx_conn)
    % connection index
    idx=anom_nodes.idx_conn(i);
    
    % signal based on cohen's D (aka snr)
    SIGNAL=snr*sqrt(ConnVarVec(idx));
    
    % add SIGNAL on the connection between the anomalous node clusters
    X_DS_tr(:,idx)=X_DS_tr(:,idx) + repmat(SIGNAL,[nDS_tr,1]);
    X_DS_ts(:,idx)=X_DS_ts(:,idx) + repmat(SIGNAL,[nDS_ts,1]);
end
% inverse fisher tx to map to (-1,+1) correlation space
X_DS_tr = tak_fisher_itransformation(X_DS_tr);
X_DS_ts = tak_fisher_itransformation(X_DS_ts);

% concatenate X_DS and X_HC to get our training and testing design matrix
Xtr=[X_DS_tr; X_HC_tr];
Xts=[X_DS_ts; X_HC_ts];

% label vector for our training data
ytr=[+ones(nDS_tr,1); -ones(nHC_tr,1)]; 
yts=[+ones(nDS_ts,1); -ones(nHC_ts,1)]; 

%==========================================================================
% pool relevant info into a struct
%==========================================================================
sim_info.nDS_tr=nDS_tr;
sim_info.nHC_tr=nHC_tr;
sim_info.ntr=ntr;

sim_info.nDS_ts=nDS_ts;
sim_info.nHC_ts=nHC_ts;
sim_info.nts=nts;

sim_info.snr=snr;
sim_info.seedPoint=seedPoint;
sim_info.anom_nodes=anom_nodes;
%% save
if fsave
    timeStamp=tak_timestamp;
    mFileName=mfilename;
    save(outPath,outVars{:})
end