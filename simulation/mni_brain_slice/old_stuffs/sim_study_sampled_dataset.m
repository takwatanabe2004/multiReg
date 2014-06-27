% sim_study_sampled_dataset (01/07/2014)
% - study the dataset sampled from sim_save_synth_fc_dataset.m
%%
clear
purge
%%
snr=0.8;
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

anom_nodes=sim_info.anom_nodes;
% load('graph_info347_2d.mat', 'coord')
load('Results_tak.mat', 'edgemaskVec','nodemask')


mean_DS_tr = mean( Xtr(ytr==+1,:));
mean_HC_tr = mean( Xtr(ytr==-1,:));
mean_DSmat_tr = tak_dvecinv(mean_DS_tr,0);
mean_HCmat_tr = tak_dvecinv(mean_HC_tr,0);

mean_DS_ts = mean( Xts(yts==+1,:));
mean_HC_ts = mean( Xts(yts==-1,:));
mean_DSmat_ts = tak_dvecinv(mean_DS_ts,0);
mean_HCmat_ts = tak_dvecinv(mean_HC_ts,0);

%% plot sample-mean-DS and sample-mean-HC figures (both from training and testing)
% figure,imexpb
% subplot(131),tplot(mean_DS_tr)
% subplot(132),tplot(mean_HC_tr)
% subplot(133),tplot(mean_DS_tr-mean_HC_tr)

figure,imexpb
subplot(131),imcov(mean_DSmat_tr)
subplot(132),imcov(mean_HCmat_tr)
subplot(133),imcov(mean_DSmat_tr-mean_HCmat_tr)
 
% figure,imexpb
% subplot(131),tplot(mean_DS_ts)
% subplot(132),tplot(mean_HC_ts)
% subplot(133),tplot(mean_DS_ts-mean_HC_ts)

figure,imexpb
subplot(131),imcov(mean_DSmat_ts)
subplot(132),imcov(mean_HCmat_ts)
subplot(133),imcov(mean_DSmat_ts-mean_HCmat_ts)
%%
anom_nodes_mask=anom_nodes.mask;
thresholded = abs(mean_DSmat_ts-mean_HCmat_ts)>0.1;
figure,imexpb
subplot(121),imedge(anom_nodes_mask)
subplot(122),imedge(thresholded)
%% plot few realizations of fc's
figure,imexpb
subplot(121),imcov( tak_dvecinv(Xtr(1,:),0))
subplot(122),imcov( tak_dvecinv(Xtr(end,:),0))

% some real-data samples
load rcorr_design_censor X
tmp1=tak_dvecinv(X(1,:),0);
tmp2=tak_dvecinv(X(2,:),0);
realdata_samp1=tmp1(nodemask,nodemask);
realdata_samp2=tmp2(nodemask,nodemask);
figure,imexpb
subplot(121),imcov(realdata_samp1)
subplot(122),imcov(realdata_samp2)
drawnow
